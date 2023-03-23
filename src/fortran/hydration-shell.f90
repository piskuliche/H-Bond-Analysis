module hyd_parameters
  implicit none
  real, parameter :: pi=DACOS(-1.D0)
  ! Main Program
  integer :: number_of_frames, number_of_atoms

  ! Read_Input
  integer :: num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  integer :: num_mol_wat, num_mol_acc
  integer :: num_other, start_wat_index, start_acc_index
  integer :: frame_start,frame_stop
  integer :: num_C, num_O, num_N, num_P, num_heavy

  real    :: dr
  real    :: critsq_C, critsq_O, critsq_P, critsq_N

  character(len=40) :: fname,iname

  ! Map_Sites
  integer :: num_lip_atoms

  ! Separator -------------
  ! Main Program
  common number_of_frames, number_of_atoms

  ! Read_Input
  common  num_mol_types, which_is_wat, which_is_acc, num_acc_sites_per_mol
  common  num_mol_wat, num_mol_acc
  common  num_other, start_wat_index, start_acc_index
  common  frame_start,frame_stop
  common  num_C, num_O, num_N, num_P, num_heavy

  common  dr
  common  critsq_C, critsq_O, critsq_P, critsq_N
  
  common  fname, iname
  ! Map_Sites
  common  num_lip_atoms
  
end module hyd_parameters

program hyd_shell

  use gmxfort_trajectory
  use hyd_parameters
  use myfuncs
  implicit none

  type (Trajectory) :: trj

  integer :: i, j, ha, wc 
  integer :: ntmp, nbins, atmp, oi

  integer :: Pocc, Cocc, Oocc, Nocc

  integer, allocatable :: atom_map(:,:)

  real :: roxsq

  real, dimension(3) :: coord
  real, dimension(4) :: criteria
  real, dimension(3,3) :: box

  real, allocatable :: rO(:,:)
  real, allocatable :: rL(:,:)
  real, allocatable :: curr_hydr(:,:)


  ! Open Log File
  open(20,file='hydration_shell.log')
  write(20,*) "Beginning gofr calculation"
  open(21,file='hydration-list.dat')
  ! Read Input File
  call Hyd_Input(criteria)
  
  ! Call Allocations
  allocate(rO(num_mol_wat,3))
  allocate(rL(num_heavy,3))
  allocate(atom_map(num_heavy,2))
  allocate(curr_hydr(num_mol_wat,2))

  ! Maps
  call Map_Sites(atom_map)

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

    Cocc = 0
    Oocc = 0
    Pocc = 0
    Nocc = 0
    !$OMP DO PRIVATE(ha,wc,roxsq) reduction(+:Pocc,Oocc,Cocc,Nocc) schedule(DYNAMIC)
    waters: do wc=1, num_mol_wat
      lips: do ha=1, num_heavy
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
      atmp = atom_map(int(curr_hydr(wc,1)),1)
      if (atmp == 1) then
        Cocc = Cocc + 1
      else if (atmp == 2) then
        Oocc = Oocc + 1
      else if (atmp == 3) then
        Pocc = Pocc + 1
      else if (atmp == 4) then
        Nocc = Nocc + 1
      endif
        
    enddo waters
    !$OMP END DO
    !$OMP SINGLE
    write(20,*) i, Cocc, Oocc, Pocc, Nocc
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

subroutine Map_Sites(atom_map)
  use hyd_parameters
  implicit none

  integer :: i, ident
  integer :: cnt, tmp

  integer, dimension(num_heavy,2) :: atom_map
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
  write(*,*) "Atoms mapped"

end subroutine Map_Sites

subroutine Hyd_Input(criteria)
  use hyd_parameters
  implicit none

  integer :: i, nother, apermol
  real    :: crit_C, crit_O, crit_P, crit_N

  real, dimension(4) :: criteria

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
  read(10,*)
  read(10,*) crit_C, crit_O, crit_P, crit_N
  criteria(1) = crit_C*crit_C
  criteria(2) = crit_O*crit_O
  criteria(3) = crit_P*crit_P
  criteria(4) = crit_N*crit_N
  

  write(*,*) "There are ", num_mol_wat, "waters"
  write(*,*) "There are ", num_mol_acc, "lipids"
  close(10)
  open(11,file='heavy.counts',status='old')
  read(11,*) num_C, num_O, num_P, num_N
  num_heavy = num_C + num_O + num_P + num_N
  close(11)
  
end subroutine Hyd_Input

