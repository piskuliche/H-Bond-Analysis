module hydr_module
  implicit none
  public
contains

  subroutine Find_Hydr(r1, r2, num_r2, comp_num, L, atom_map, criteria, hydration)
    ! Finds the closest water to atoms in the given group
    
    ! Args:
    ! r1 - coordinates of the atoms in the first group (3)
    ! r2 - coordinates of the atoms in the second group (num_r2,3)
    ! num_r2 - number of atoms in the second group
    ! comp_num - the component number of the second group
    ! L - box dimensions (3,3)
    ! atom_map - mapping of the atoms in the second group to the first (num_r2,2)
    ! criteria - criteria for the atoms in the second group (num_types - assumed)
    ! hydration - the atom, and square distance of closest atom (2)

    ! Returns:
    ! hydration - the atom, and square distance of closest atom (2)

    ! Local Variables
    ! i - loop variable
    ! roxsq - squared distance between the two atoms



    use myfuncs
    implicit none
    integer, intent(in) :: num_r2, comp_num
    real, dimension(:), intent(in) :: r1
    real, dimension(:,:), intent(in) :: r2, L
    integer, dimension(:,:), intent(in) :: atom_map
    real, dimension(:,:), intent(in) :: criteria
    real, dimension(:), intent(inout) :: hydration
 
    integer :: i
    real :: roxsq 
    ! Loop over all the atoms in the second group
    do i = 1, num_r2
      ! Calculate the squared distance between the two atoms (with PBC)
      roxsq = distance2(r1, r2(i,:), L)
      if (roxsq <= criteria(comp_num, atom_map(i,2))) then
        if (hydration(2) > 0.0) then
          if (roxsq < hydration(2)) then
            hydration(1) = i
            hydration(2) = roxsq
          endif
        else
          hydration(1) = i
          hydration(2) = roxsq
        endif 
      endif
    end do

  end subroutine Find_Hydr


  subroutine Alt_Hyd_Input(mapfile, frame_start, frame_stop, fname, iname &
                            , num_components, num_heavy, component_start &
                            , atoms_per_component &
                            , is_water, criteria)
    implicit none

    integer :: i, j
    integer :: frame_start, frame_stop
    integer :: num_components, atom_count
    integer :: is_water
    integer :: check_if_water

    integer, allocatable :: component_label(:), num_mol(:), atoms_per_mol(:), num_acc_sites_per_mol(:)
    integer, allocatable :: component_start(:), is_mem_component(:), num_heavy(:)
    integer, allocatable :: atoms_per_component(:)
    
    real, allocatable :: criteria(:,:), crit(:,:)

    character(len=40) :: fname, iname, mapfile

    write(*,*) "Reading in the input file"

    open(10,file='hydration.in',status='old')
    read(10,*)
    read(10,*) fname, iname, mapfile
    read(10,*)
    read(10,*) frame_start,frame_stop
    read(10,*)
    read(10,*) num_components
    read(10,*)

    ! Allocate Component based arrays
    allocate(component_label(num_components))
    allocate(num_mol(num_components))
    allocate(atoms_per_mol(num_components))
    allocate(num_acc_sites_per_mol(num_components))
    allocate(component_start(num_components))
    allocate(is_mem_component(num_components))
    allocate(crit(num_components,4))
    allocate(criteria(num_components,4))
    allocate(num_heavy(num_components))
    allocate(atoms_per_component(num_components))

    component_label=0; num_mol=0; atoms_per_mol=0; num_acc_sites_per_mol=0
    component_start=0; is_mem_component=0; criteria=0; num_heavy=0
    atoms_per_component=0

    ! This reads in the component labels, number of molecules, number of acceptor sites,
    ! and whether or not the component is in membrane or is a water molecule
    ! also calculates start of the components
    atom_count = 0
    is_water=-1
    do i=1, num_components
        check_if_water = -1
        ! Read the Component Information
        read(10,*) component_label(i), num_mol(i), atoms_per_mol(i), num_acc_sites_per_mol(i)
        read(10,*) check_if_water

        ! Check if the component is a water molecule
        if ( check_if_water == 1 ) then
            ! Error test for if there is more than one water component
            if ( is_water /= -1 ) then
                write(*,*) "Error: More than one water component found!"
                write(*,*) "To solve this error, run the code individually"
                write(*,*) "For each water component"
            endif
            if ( num_acc_sites_per_mol(i) /= 1 ) then
                write(*,*) "Error: Water molecules must have one site per molecule!"
                write(*,*) "This should be placed on the Oxygen atom"
                stop
            endif
            is_water = i
        endif

        write(*,*) "There are ", num_mol(i), " ", component_label(i), " molecules"
        ! Do bookkeeping for the number of atoms and start of each component.
        component_start(i) = atom_count
        ! Number of heavy atoms is the number of acceptor sites per molecule times the number of molecules
        num_heavy(i) = num_acc_sites_per_mol(i)*num_mol(i)
        ! Number of atoms per component is the number of molecules times the number of atoms per molecule
        atoms_per_component(i) = num_mol(i)*atoms_per_mol(i)
        ! Total number of atoms is the sum of the atoms per component
        atom_count = atom_count + atoms_per_component(i)
    enddo

    read(10,*)
    

    ! Read in the solvation shell criteria
    crit = 0.0; criteria = 0.0
    do i=1, num_components
      read(10,*) crit(i,1), crit(i,2), crit(i,3), crit(i,4)
    enddo
    
    ! Square the criteria and store
    do i=1, num_components
      do j=1, 4
        criteria(i,j) = crit(i,j)*crit(i,j)
      enddo
    enddo
    
    close(10)
    
    deallocate(crit)

    write(*,*) "Finished reading in the input file"
    
  end subroutine Alt_Hyd_Input

end module hydr_module