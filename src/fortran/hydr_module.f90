module hydr_module
  implicit none
  public
contains

  subroutine Find_Hydr(r1, r2, num_r2, L, atom_map, criteria, hydration)
    ! Finds the closest water to atoms in the given group
    
    ! Args:
    ! r1 - coordinates of the atoms in the first group (3)
    ! r2 - coordinates of the atoms in the second group (num_r2,3)
    ! num_r2 - number of atoms in the second group
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
    integer, intent(in) :: num_r2
    real, dimension(:), intent(in) :: r1
    real, dimension(:,:), intent(in) :: r2, L
    integer, dimension(:,:), intent(in) :: atom_map
    real, dimension(:), intent(in) :: criteria
    real, dimension(:), intent(inout) :: hydration
 
    integer :: i
    real :: roxsq

    ! Loop over all the atoms in the second group
    do i = 1, num_r2
      ! Calculate the squared distance between the two atoms (with PBC)
      roxsq = distance2(r1, r2(i,:), L)
      if (roxsq <= criteria(atom_map(i,2))) then
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

  subroutine Generate_Atom_Map(num_atoms, component_start, num_heavy, atom_map)
    ! This subroutine maps the lipid/membrane atoms that are heavy and their type.
    ! This is used to map all the coordinates in the system in terms of only the heavy atoms
    ! Args:
    !   num_atoms - number of atoms in each component
    !   component_start - the starting index of each component
    !   num_heavy - the number of heavy atoms in each component
    !   atom_map - the mapping of the heavy atoms to the atom index and type
    ! Returns:
    !   atom_map - the mapping of the heavy atoms to the atom index and type

    ! Input Files:
    !   map_of_heavy_atoms.info - the mapping of the heavy atoms to the atom index and type
    !   Note - that right now, four atom types of interest are assumed for each component
    !   Which are based on the criteria array (see Alt_Hyd_Input) 
    !   For each component, label 0 for any that isn't heavy (e.g. hydrogen)
    !   and (1-4) + j*4 for the heavy atoms of interest, where j is the component number
    !   Note - that num_heavy atoms must be equivalent to the number of acceptor sites per molecule
    !   Times the number of molecules that make up that component.


    implicit none
    integer :: i, j, cnt
    integer :: atom_type
    integer, dimension(:), intent(in) :: num_heavy, num_atoms, component_start
    integer, dimension(:,:,:), intent(inout) :: atom_map

    atom_map = 0
    open(12, file='map_of_heavy_atoms.info', status='old')
    do i = 1, size(num_heavy)
        cnt = 0
        do j = 1, num_atoms(i)
            read(12,*) atom_type
            if (atom_type /= 0) then
                atom_map(i,cnt,1) = j + component_start(j)
                atom_map(i,cnt,2) = atom_type
                cnt = cnt + 1
                endif
        enddo
        if ( cnt /= num_heavy(i) ) then
            write(*,*) "Error: Number of heavy atoms does not match the number of acceptor sites per molecule"
            write(*,*) "Number of heavy atoms: ", cnt
            write(*,*) "Number of acceptor sites per molecule: ", num_heavy(i)
            write(*,*) "Component: ", i
            stop
        endif
    enddo 

  end subroutine Generate_Atom_Map


  subroutine Alt_Hyd_Input(frame_start, frame_stop, fname, iname &
                            , num_components, num_heavy, component_start &
                            , atoms_per_component &
                            , is_water, criteria)
    implicit none

    integer :: i
    integer :: frame_start, frame_stop
    integer :: num_components, atom_count
    integer :: is_water
    integer :: check_if_water


    real :: dr

    integer, allocatable :: component_label(:), num_mol(:), atoms_per_mol(:), num_acc_sites_per_mol(:)
    integer, allocatable :: component_start(:), is_mem_component(:), num_heavy(:)
    integer, allocatable :: atoms_per_component(:)
    
    real, allocatable :: criteria(:), crit(:)

    character(len=40) :: fname, iname

    open(10,file='hydration.in',status='old')
    read(10,*)
    read(10,*) fname, iname
    read(10,*)
    read(10,*) frame_start,frame_stop
    read(10,*)
    read(10,*) dr, num_components
    read(10,*)

    ! Allocate Component based arrays
    allocate(component_label(num_components))
    allocate(num_mol(num_components))
    allocate(atoms_per_mol(num_components))
    allocate(num_acc_sites_per_mol(num_components))
    allocate(component_start(num_components))
    allocate(is_mem_component(num_components))
    allocate(crit(4*num_components))
    allocate(criteria(4*num_components))
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
      read(10,*) crit((i-1)*4+1), crit((i-1)*4+2), crit((i-1)*4+3), crit((i-1)*4+4)
    enddo
    
    ! Square the criteria and store
    do i=1, 4*num_components
      criteria(i) = crit(i)*crit(i)
    enddo
    
    close(10)
    
    deallocate(crit)
    
  end subroutine Alt_Hyd_Input

end module hydr_module