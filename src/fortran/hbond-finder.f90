module constants
  implicit none
  real,parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: number_of_frames, number_of_atoms


  ! Read_Input
  integer :: num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  integer :: num_mol_wat, num_mol_acc, num_acc_sites
  integer :: num_other, start_wat_index

  real :: thickness
  real :: crit_ox, crit_hx, crit_ohx
  real :: critsq_ox, critsq_hx

  character(len=40) :: fname,iname

  ! Map Acceptor
  integer :: num_acc_atoms

  ! Separator -----------------------------

  common  number_of_frames, number_of_atoms

  ! Read_Input
  common  num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  common  num_mol_wat, num_mol_acc, num_acc_sites
  common  num_other, start_wat_index

  common  thickness
  common  crit_ox, crit_hx, crit_ohx
  common  critsq_ox, critsq_hx

  common  fname, iname

  ! Map Acceptor
  common num_acc_atoms

end module parameters


program hba
  use gmxfort_trajectory
  use constants
  use parameters
  use myfuncs
  implicit none

  type (Trajectory) :: trj

  integer :: i, j, ac, wc,oh,w2
  integer :: ntmp,oi,hi,found

  real :: roxsq, rhxsq, thta
  real :: nhbonds,hb2

  integer, allocatable :: acc_map(:)
  integer, allocatable :: prev_hbond(:), curr_hbond(:)
  
  real, dimension(3) :: coord
  real, dimension(3,3) :: box

  real, allocatable :: r_water(:,:,:)
  real, allocatable :: r_acceptor(:,:)



  ! Read Input File
  call Read_Input()

  ! Call allocations
  allocate(r_water(num_mol_wat,3,3))
  allocate(r_acceptor(num_acc_sites,3))
  allocate(acc_map(num_acc_sites))
  allocate(curr_hbond(num_mol_wat*2))
  allocate(prev_hbond(num_mol_wat*2))

  ! Generate Map of Acceptor Sites
  call Map_Acceptors(acc_map)
  start_wat_index = num_acc_atoms+num_other

  ! Open Trajectory
  call trj%open(trim(fname),trim(iname))
  number_of_frames = trj%nframes
  number_of_atoms  = trj%natoms()
  write(*,*) "There are ", number_of_frames," frames."
  write(*,*) "There are ", number_of_atoms," atoms."

  ! Loop over frames
  do i=1,number_of_frames
    ! Read Frame
    ntmp = trj%read_next(1)
    box  = trj%box(1)
    ! Assign Acc Atoms
    do ac=1, num_acc_sites
      coord(:) = trj%x(1,acc_map(ac))
      r_acceptor(ac,:) = coord(:)
    enddo
    ! Assign Water Atoms
    write(*,*) start_wat_index
    do wc=1,num_mol_wat
      oi = start_wat_index + (wc-1)*3
      do j=1,3
        coord(:) = trj%x(1,oi+j)
        r_water(wc,j,:) = coord(:)
      enddo
    enddo
    
    ! Call Distance Calculation
    nhbonds = 0
    hb2=0
    prev_hbond=0
    curr_hbond=0
    !$OMP PARALLEL DO PRIVATE(oh,ac,w2,found,roxsq,rhxsq,thta) reduction(+:nhbonds,hb2)
    waters: do wc=1,num_mol_wat
      ohs: do oh=1,2
        found=0
        accs: do ac=1,num_acc_sites
          roxsq = distance2(r_water(wc,1,:),r_acceptor(ac,:),box)
          if (roxsq <= critsq_ox) then
            rhxsq = distance2(r_water(wc,oh+1,:),r_acceptor(ac,:),box)
            if (rhxsq <= critsq_hx) then
              thta = bond_angle(r_water(wc,oh+1,:),r_water(wc,1,:),r_acceptor(ac,:),box)
              if (thta <= crit_ohx) then
                nhbonds = nhbonds + 1
                curr_hbond((wc-1)*2+oh)=ac
                found=1
                exit accs
              endif
            endif
          endif
        enddo accs
        if (found == 0) then
          wat2: do w2=1,num_mol_wat
            roxsq = distance2(r_water(wc,1,:),r_water(w2,1,:),box)
            if (roxsq <= critsq_ox) then
              rhxsq = distance2(r_water(wc,oh+1,:),r_water(w2,1,:),box)
              if (rhxsq <= critsq_hx) then
                thta = bond_angle(r_water(wc,oh+1,:),r_water(wc,1,:),r_water(w2,1,:),box)
                if (thta <= crit_ohx) then
                  hb2 = hb2 + 1
                  curr_hbond((wc-1)*2+oh)=w2+num_acc_sites
                  exit wat2
                endif !thta
              endif ! rhxsq
            endif ! Roxsq
          enddo wat2
        endif ! found
      enddo ohs
    enddo waters
    !$OMP END PARALLEL DO
  
    write(*,*) nhbonds,hb2
  ! End frame loop
  enddo
  
  ! End Program

end program hba

subroutine Read_Input()
  use parameters
  use constants
  implicit none

  integer :: i
  integer :: nother

  character(len=6) :: ctmp


  nother = 0
  open(10,file='hbonds.in',status='old')
  read(10,*)
  read(10,*) fname,iname
  read(10,*)
  read(10,*) thickness
  read(10,*)
  read(10,*) crit_ox, crit_hx, crit_ohx
  critsq_ox = crit_ox*crit_ox
  critsq_hx = crit_hx*crit_hx
  read(10,*)
  read(10,*) num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  read(10,*)
  do i=1,num_mol_types
    if (i == which_is_wat) then
      read(10,*) ctmp, num_mol_wat
    else if (i == which_is_acc) then
      read(10,*) ctmp, num_mol_acc
    else
      read(10,*) ctmp, nother 
      num_other = num_other + nother
    endif
  enddo
  write(*,*) "There are ", num_mol_wat, " waters"
  ! Calculate number of acceptor sites
  num_acc_sites = num_acc_sites_per_mol * num_mol_acc
  close(10)
end subroutine

subroutine Map_Acceptors(map)
  use parameters
  use constants
  implicit none

  integer :: i
  
  integer :: cnt,tmp

  integer, dimension(num_acc_sites) :: map

  cnt = 1
  
  open(11,file='is_acc.txt',status='old')
  read(11,*) tmp, num_acc_atoms
  do i=1,num_acc_atoms
    read(11,*) tmp
    if (tmp == 1) then
      map(cnt) = i
      cnt = cnt +  1
    endif
  enddo
  close(11)
end subroutine Map_Acceptors

