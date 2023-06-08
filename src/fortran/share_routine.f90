module share_routine
  implicit none
  public
contains
  subroutine Generate_Atom_Map(filename, num_atoms, component_start, num_heavy, atom_map)
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

      character(len=40), intent(in) :: filename

      write(*,*) "Generating the atom map"

      atom_map = 0
      open(12, file=trim(filename), status='old')
      do i = 1, size(num_heavy)
          write(*,*) "Mapping Component: ", i
          write(*,*) "With ", num_heavy(i), " heavy atoms"
          write(*,*) "And ", num_atoms(i), " total atoms"
          write(*,*) "Component starts at atom: ", component_start(i)
          cnt = 0
          do j = 1, num_atoms(i)
            if ( cnt > num_heavy(i) ) then
              write(*,*) "Error: Number of heavy atoms does not match the number of acceptor sites per molecule"
              write(*,*) "Number of heavy atoms: ", cnt
              write(*,*) "Number of acceptor sites per molecule: ", num_heavy(i)
              write(*,*) "Component: ", i
              stop
            endif
            read(12,*) atom_type
            if (atom_type /= 0) then
              cnt = cnt + 1
              atom_map(i,cnt,1) = j + component_start(i)
              atom_map(i,cnt,2) = atom_type
            endif
          enddo
          
      enddo 

      write(*,*) "Finished generating the atom map"

  end subroutine Generate_Atom_Map
end module share_routine
