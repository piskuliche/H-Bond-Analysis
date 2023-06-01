module hyd_parameters
  implicit none
  real, parameter :: pi=DACOS(-1.D0)
  ! Main Program
  integer :: number_of_frames, number_of_atoms

  ! Read_Input
  integer :: num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  integer :: num_mol_wat, num_mol_acc
  integer :: num_other, start_wat_index, start_acc_index, start_laur_index
  integer :: num_mol_laur
  integer :: num_acc_sites_per_laur, which_is_laur
  integer :: frame_start,frame_stop
  integer :: num_heavy, num_laur_heavy

  real    :: dr
  real    :: critsq_C, critsq_O, critsq_P, critsq_N

  character(len=40) :: fname,iname

  ! Map_Sites
  integer :: num_lip_atoms, num_laur_atoms

  ! Separator -------------
  ! Main Program
  common number_of_frames, number_of_atoms

  ! Read_Input
  common  num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  common  num_mol_wat, num_mol_acc
  common  num_other, start_wat_index, start_acc_index, start_laur_index
  common  num_mol_laur
  common  num_acc_sites_per_laur, which_is_laur
  common  frame_start,frame_stop
  common  num_heavy, num_laur_heavy

  common  dr
  common  critsq_C, critsq_O, critsq_P, critsq_N
  
  common  fname, iname
  ! Map_Sites
  common  num_lip_atoms, num_laur_atoms
  
end module hyd_parameters

program hyd_shell

  use gmxfort_trajectory
  use hyd_parameters
  use myfuncs
  implicit none

  type (Trajectory) :: trj

  integer :: i, j, ha, wc 
  integer :: ntmp, nbins, atmp, oi
  integer :: Pocc, Oocc, Cocc, Nocc
  integer :: Pocc2, Oocc2, Cocc2, Nocc2
  integer :: num_total_heavy

  integer, allocatable :: atom_map(:,:), laur_atom_map(:,:)

  real :: roxsq

  real, dimension(3) :: coord
  real, dimension(8) :: criteria
  real, dimension(3,3) :: box

  real, allocatable :: rO(:,:)
  real, allocatable :: rL(:,:)
  real, allocatable :: curr_hydr(:,:)


  ! Open Log File
  open(20,file='hydration_shell.log')
  write(20,*) "Beginning hydration calculation"
  open(21,file='hydration-list.dat')
  
  ! Read Input File
  call Hyd_Input(criteria)

  write(*,*) "Criteria are: ", criteria
  
  ! Call Allocations
  allocate(rO(num_mol_wat,3))
  allocate(rL(num_heavy+num_laur_heavy,3))
  allocate(atom_map(num_heavy+num_laur_heavy,2))
  allocate(curr_hydr(num_mol_wat,2))

  ! Maps
  call Map_Acc_Sites(atom_map)

  ! Open Trajectory
  call trj%open(trim(fname),trim(iname))
  number_of_frames = trj%nframes
  number_of_atoms  = trj%natoms()


  write(20,*) "There are ", number_of_frames," frames."
  write(20,*) "There are ", number_of_atoms," atoms."

  do i=1, frame_start-1
    ntmp = trj%read_next(1)
    if (mod(i,100) == 0) write(*,*) "Cycled through frame ", i
  enddo
  !$OMP PARALLEL
  do i=frame_start,frame_stop
    !$OMP SINGLE
    write(*,*) "Reached frame ", i
    ! Read Frame
    ntmp = trj%read_next(1)
    box  = trj%box(1)
    ! Assign Lip Atoms
    do ha=1,num_heavy
      coord(:) = trj%x(1,atom_map(ha,1))
      rL(ha,:) = coord(:)
    enddo
    if (num_laur_heavy > 0) then
      do ha=num_heavy, num_heavy+num_laur_heavy
        coord(:) = trj%x(1,atom_map(ha,1))
        rL(ha,:) = coord(:)
      enddo
    endif
    ! Assign Water Atoms
    do wc=1,num_mol_wat
      oi = start_wat_index + (wc-1)*3
      coord(:) = trj%x(1,oi+1)
      rO(wc,:) = coord(:)
    enddo


    ! Distance Calculation
    curr_hydr = 0.0
    !$OMP END SINGLE
    !$OMP BARRIER

    Cocc=0; Oocc=0; Pocc=0; Nocc=0
    Cocc2=0; Oocc2=0; Pocc2=0; Nocc2=0
    num_total_heavy = num_heavy + num_laur_heavy
    !$OMP DO PRIVATE(ha,wc,roxsq) reduction(+:Pocc,Oocc,Cocc,Nocc) schedule(DYNAMIC)
    waters: do wc=1, num_mol_wat
      lips: do ha=1, num_total_heavy
        roxsq = distance2(rO(wc,:),rL(ha,:),box)
        ! Check that it meets criteria
        if (roxsq <= criteria(atom_map(ha,2))) then
          ! Check that ther is not already something found
          if (curr_hydr(wc,2) > 0.0) then
            ! Check distance is smallest, if it is then save
            if (roxsq < curr_hydr(wc,2)) then
              curr_hydr(wc,1) = ha
              curr_hydr(wc,2) = roxsq
            endif
          else ! If zero, save
            curr_hydr(wc,1) = ha
            curr_hydr(wc,2) = roxsq
          endif 
        endif
      enddo lips
      ! Map info for current water
      atmp = atom_map(int(curr_hydr(wc,1)),2)
      if (atmp == 1) then
        Cocc = Cocc + 1
      else if (atmp == 2) then
        Oocc = Oocc + 1
      else if (atmp == 3) then
        Pocc = Pocc + 1
      else if (atmp == 4) then
        Nocc = Nocc + 1
      else if (atmp == 5) then
        Cocc2 = Cocc2 + 1
      else if (atmp == 6) then
        Oocc2 = Oocc2 + 1
      else if (atmp == 7) then
        Pocc2 = Pocc2 + 1
      else if (atmp == 8) then
        Nocc2 = Nocc2 + 1
      endif
        
    enddo waters
    !$OMP END DO
    !$OMP SINGLE
    write(20,*) i, Cocc, Oocc, Pocc, Nocc, Cocc2, Oocc2, Pocc2, Nocc2
    do wc=1,num_mol_wat
      write(21,*) int(curr_hydr(wc,1)), curr_hydr(wc,2)
      write(21,*) int(curr_hydr(wc,1)), curr_hydr(wc,2)
    enddo 
    !$OMP END SINGLE
  enddo !frame
  !$OMP END PARALLEL
  write(20,*) "Calculation complete"
  close(20)
  close(21)


end program hyd_shell

subroutine Map_Acc_Sites(atom_map)
  use hyd_parameters
  implicit none

  integer :: i, ident
  integer :: cnt, tmp

  integer, dimension(num_heavy+num_laur_heavy,2) :: atom_map
  write(*,*) num_heavy
  write(*,*) "Mapping atoms"
  open(12,file='is_heavy.txt',status='old')
  read(12,*) num_lip_atoms
  cnt = 1
  do i=1,num_lip_atoms
    read(12,*) tmp
    if (tmp /= 0) then
      atom_map(cnt,1) = i + start_acc_index
      atom_map(cnt,2) = tmp
      cnt = cnt + 1
    endif
  enddo
  close(12)
  ! Repeat for laurdan if needed
  if (num_mol_laur > 0) then
    open(12,file='laur_is_heavy.txt',status='old')
    read(12,*) num_laur_atoms
    do i=1,num_laur_atoms
      read(12,*) tmp
      if (tmp /= 0) then
        atom_map(cnt,1) = i + start_laur_index
        atom_map(cnt,2) = tmp
        cnt = cnt + 1
      endif
    enddo
    close(12)
  endif
  write(*,*) "Atoms mapped"

end subroutine Map_Acc_Sites



subroutine Hyd_Input(criteria)
  use hyd_parameters
  implicit none

  integer :: i, nother, apermol
  integer :: atom_count
  integer :: nmols
  integer :: num_C, num_P, num_O, num_N
  real    :: crit_C, crit_O, crit_P, crit_N

  real, dimension(8) :: criteria

  character(len=6) :: ctmp

  open(10,file='hbonds.in',status='old')
  read(10,*)
  read(10,*) fname, iname
  read(10,*)
  read(10,*) frame_start,frame_stop
  read(10,*)
  read(10,*) dr
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*) num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  read(10,*)
  read(10,*) which_is_laur, num_acc_sites_per_laur
  read(10,*)
  ! Loop for reading in number of molecules and atoms per molecule
  ! Also calculate the starting index for the water molecules
  ! and the starting index for the acceptor molecules
  start_wat_index = 0
  start_acc_index = 0
  start_laur_index= 0
  num_mol_acc = 0; num_mol_laur = 0; num_mol_wat = 0
  atom_count=0
  do i=1, num_mol_types
    read(10,*) ctmp, nmols, apermol
    if (i == which_is_wat) then
      num_mol_wat = nmols
      start_wat_index = atom_count
    endif
    if (i == which_is_acc) then
      num_mol_acc = nmols
      start_acc_index = atom_count
    endif
    if (i == which_is_laur) then
      num_mol_laur = nmols
      start_laur_index = atom_count
    endif
    atom_count = atom_count + nmols*apermol
  enddo

  read(10,*)
  read(10,*) crit_C, crit_O, crit_P, crit_N
  criteria(1) = crit_C*crit_C
  criteria(2) = crit_O*crit_O
  criteria(3) = crit_P*crit_P
  criteria(4) = crit_N*crit_N
  do i=5,8
    criteria(i) = criteria(i-4)
  enddo
  

  write(*,*) "There are ", num_mol_wat, "waters"
  write(*,*) "There are ", num_mol_acc, "lipids"
  write(*,*) "There are ", num_mol_laur, "laurdans"
  close(10)
  ! Read in heavy atoms
  open(11,file='heavy.counts',status='old')
  num_C = 0; num_O = 0; num_P = 0; num_N = 0
  read(11,*) num_C, num_O, num_P, num_N
  num_heavy = num_C + num_O + num_P + num_N
  close(11)
  num_C = 0; num_O = 0; num_P = 0; num_N = 0
  open(11,file='laur.counts',status='old')
  read(11,*) num_C, num_O, num_P, num_N
  num_laur_heavy = num_C + num_O + num_P + num_N
  close(11)
  
end subroutine Hyd_Input

