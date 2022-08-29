module share_routine
  implicit none
  public
contains
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
    read(10,*) frame_start, frame_stop
    read(10,*)
    read(10,*) thickness
    read(10,*)
    read(10,*) crit_ox, crit_hx, crit_deg
    critsq_ox = crit_ox*crit_ox
    critsq_hx = crit_hx*crit_hx
    crit_rad = crit_deg*pi/180.0
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
end module share_routine
