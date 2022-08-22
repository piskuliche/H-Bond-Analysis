module constants
  implicit none
  real,parameter :: pi=DACOS(-1.D0)
end module constants

module parameters
  integer :: number_of_frames, number_of_atoms

  ! Read_Input
  real :: thickness
  real :: crit_ox, crit_hx, crit_ohx


  character(len=40) :: fname,iname


  common :: number_of_frames, number_of_atoms

  ! Read_Input
  common :: thickness
  common :: crit_ox, crit_hx, crit_ohx

  common :: fname, iname

end module parameters

program hba
  use gmxfort_trajectory
  use constants
  use parameters
  implicit none

  type (Trajectory) :: trj

  integer :: i

  ! Read Input File
  call Read_Input()

  ! Open Trajectory
  call trj%open(trim(fname),trim(iname))
  number_of_frames = trj%nframes
  number_of_aoms  = trj%natoms()
  write(*,*) "There are ", number_of_frames," frames."
  write(*,*) "There are ", number_of_atoms," atoms."

  ! Loop over frames
  do i=1,number_of_frames
    call Find_Hbonds()
  enddo
  
  ! End Program

end program hba

subroutine Read_Input()
  use parameters
  use constants
  open(10,file='hbonds.in',status='old')
  read(10,*)
  read(10,*) fname,iname
  read(10,*)
  read(10,*) thickness
  read(10,*)
  read(10,*) crit_ox, crit_hx, crit_ohx
  read(10,*)
  close(10)
end subroutine

subroutine Find_Hbonds()
  use parameters
  use constants
  implicit none
  write(*,*) "Hbonds"
end subroutine Find_Hbonds

