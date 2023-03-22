module constants
  implicit none
  real,parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: number_of_frames, number_of_atoms, num_ohs


  ! Read_Input
  integer :: num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  integer :: num_mol_wat, num_mol_acc, num_acc_sites
  integer :: num_other, start_wat_index, start_acc_index
  integer :: frame_start,frame_stop

  real :: thickness
  real :: crit_ox, crit_hx, crit_deg
  real :: critsq_ox, critsq_hx, crit_rad


  character(len=40) :: fname,iname

  ! Map Acceptor
  integer :: num_acc_atoms

  ! Separator -----------------------------

  common  number_of_frames, number_of_atoms, num_ohs

  ! Read_Input
  common  num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  common  num_mol_wat, num_mol_acc, num_acc_sites
  common  num_other, start_wat_index, start_acc_index
  common  frame_start,frame_stop

  common  thickness
  common  crit_ox, crit_hx, crit_deg
  common  critsq_ox, critsq_hx, crit_rad

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
  real :: acc_hbonds, wat_hbonds

  integer, allocatable :: acc_map(:)
  integer, allocatable :: prev_hbond(:), curr_hbond(:)
  
  real, dimension(3) :: coord
  real, dimension(3,3) :: box

  real, allocatable :: r_water(:,:,:)
  real, allocatable :: r_acceptor(:,:)



  ! Open Log File
  open(12,file='hbonds.log')
  write(12,*) "Beginning Hbond Calculation"
  open(13,file='hbond-list.dat')
  ! Read Input File
  call Read_Input()
  num_ohs = num_mol_wat*2

  ! Call Error Testing
  call Error_Testing()

  ! Call allocations
  allocate(r_water(num_mol_wat,3,3))
  allocate(r_acceptor(num_acc_sites,3))
  allocate(acc_map(num_acc_sites))
  allocate(curr_hbond(num_ohs))
  allocate(prev_hbond(num_ohs))

  ! Generate Map of Acceptor Sites
  call Map_Acceptors(acc_map)

  ! Open Trajectory
  call trj%open(trim(fname),trim(iname))
  number_of_frames = trj%nframes
  number_of_atoms  = trj%natoms()
  write(12,*) "There are ", number_of_frames," frames."
  write(12,*) "There are ", number_of_atoms," atoms."

  write(*,*) "Begin loop over frames"
  write(12,*) "#frame wat_hbonds acc_hbonds"
  do i=1,frame_start-1
    ntmp = trj%read_next(1)
    if (mod(i,100) == 0) write(*,*) "Cycled through frame ", i
  enddo
  ! Loop over frames
  !$OMP PARALLEL
  do i=frame_start,frame_stop
    !$OMP SINGLE
    write(*,*) "Reached frame ",i
    ! Read Frame
    ntmp = trj%read_next(1)
    box  = trj%box(1)
    ! Assign Acc Atoms
    do ac=1, num_acc_sites
      coord(:) = trj%x(1,acc_map(ac))
      r_acceptor(ac,:) = coord(:)
    enddo
    ! Assign Water Atoms
    do wc=1,num_mol_wat
      oi = start_wat_index + (wc-1)*3
      do j=1,3
        coord(:) = trj%x(1,oi+j)
        r_water(wc,j,:) = coord(:)
      enddo
    enddo
    
    ! Call Distance Calculation
    acc_hbonds = 0
    wat_hbonds = 0
    prev_hbond = 0 
    curr_hbond = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    !$OMP DO PRIVATE(oh,ac,w2,found,roxsq,rhxsq,thta) reduction(+:acc_hbonds,wat_hbonds) SCHEDULE(DYNAMIC)
    waters: do wc=1,num_mol_wat
      ohs: do oh=1,2
        found=0
        accs: do ac=1,num_acc_sites
          roxsq = distance2(r_water(wc,1,:),r_acceptor(ac,:),box)
          if (roxsq <= critsq_ox) then
            rhxsq = distance2(r_water(wc,oh+1,:),r_acceptor(ac,:),box)
            if (rhxsq <= critsq_hx) then
              thta = bond_angle(r_water(wc,oh+1,:),r_water(wc,1,:),r_acceptor(ac,:),box)
              if (thta <= crit_rad) then
                acc_hbonds = acc_hbonds + 1
                curr_hbond((wc-1)*2+oh)=ac
                found=1
                exit accs
              endif
            endif
          endif
        enddo accs
        if (found == 0) then
          wat2: do w2=1,num_mol_wat
            if ( wc /= w2) then
              roxsq = distance2(r_water(wc,1,:),r_water(w2,1,:),box)
              if (roxsq <= critsq_ox) then
                rhxsq = distance2(r_water(wc,oh+1,:),r_water(w2,1,:),box)
                if (rhxsq <= critsq_hx) then
                  thta = bond_angle(r_water(wc,oh+1,:),r_water(wc,1,:),r_water(w2,1,:),box)
                  if (thta <= crit_rad) then
                    wat_hbonds = wat_hbonds + 1
                    curr_hbond((wc-1)*2+oh)=w2+num_acc_sites
                    exit wat2
                  endif !thta
                endif ! rhxsq
              endif ! Roxsq
            endif ! check equiv
          enddo wat2
        endif ! found
      enddo ohs
    enddo waters
    !$OMP END  DO
    !$OMP SINGLE
    write(12,*) i, wat_hbonds, acc_hbonds 
    do wc=1,num_ohs
      write(13,*) curr_hbond(wc)
    enddo
    !$OMP END SINGLE
    !$OMP BARRIER
  ! End frame loop
  enddo
  !$OMP END PARALLEL
  write(12,*) "Calculation complete"
  close(12)
  
  ! End Program

end program hba

subroutine Error_Testing()
  use parameters 
  use constants
  USE iso_fortran_env, ONLY : error_unit
  implicit none
  
  ! Test that criteria are not zero
  if (crit_ox*crit_hx*crit_deg == 0) then
    write(error_unit,*) "Error: Criteria must be nonzero"
    stop 1
  endif
  if (frame_start >= frame_stop) then
    write(error_unit,*) "Error: frame_start must be less than frame_stop"
    stop 1
  endif
end subroutine Error_Testing


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
      map(cnt) = i + start_acc_index ! This handles the offset if there are other atoms before
      cnt = cnt +  1
    endif
  enddo
  close(11)
end subroutine Map_Acceptors

subroutine Read_Input()
    use parameters
    use constants
    implicit none

    integer :: i
    integer :: nother, apermol

    character(len=6) :: ctmp


    nother = 0
    apermol = 0
    ! Open and read the input file
    open(10,file='hbonds.in',status='old')
    read(10,*)
    read(10,*) fname,iname
    read(10,*)
    read(10,*) frame_start, frame_stop
    read(10,*)
    read(10,*) thickness
    read(10,*)
    read(10,*) crit_ox, crit_hx, crit_deg
    ! Do some calculations to get square cutoffs and radians
    critsq_ox = crit_ox*crit_ox
    critsq_hx = crit_hx*crit_hx
    crit_rad = crit_deg*pi/180.0
    read(10,*)
    read(10,*) num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
    read(10,*)
    ! Read in the number of each molecule type
    start_wat_index = 0
    start_acc_index = 0
    do i=1,num_mol_types
      if (i == which_is_wat) then
        read(10,*) ctmp, num_mol_wat, apermol
      else if (i == which_is_acc) then
        read(10,*) ctmp, num_mol_acc, apermol
        if (which_is_acc < which_is_wat) then !add if acc is before wat
          start_wat_index = start_wat_index + num_mol_acc*apermol
        endif
      else if (i /= which_is_wat .and. i /= which_is_acc) then
        read(10,*) ctmp, nother, apermol
        num_other = num_other + nother*apermol
        if (i < which_is_wat) then ! add if other is before wat
          start_wat_index = start_wat_index + nother*apermol
        endif
        if (i < which_is_acc) then ! add if other is before acc
          start_acc_index = start_acc_index + nother*apermol
        endif
      endif
    enddo
    write(*,*) "There are ", num_mol_wat, " waters"
    write(*,*) "There are ", num_mol_acc, " acceptors"
    write(*,*) "There are ", num_other, " other atoms"
    write(*,*) "The start index of the water molecules is ", start_wat_index
    write(*,*) "The start index of the acceptor molecules is ", start_acc_index
    ! Calculate number of acceptor sites
    num_acc_sites = num_acc_sites_per_mol * num_mol_acc
    close(10)
end subroutine Read_Input
