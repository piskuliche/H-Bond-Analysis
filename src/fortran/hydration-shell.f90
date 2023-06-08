program hyd_shell
  use hydr_module
  use gmxfort_trajectory
  use myfuncs
  use share_routine
  implicit none

  type (Trajectory) :: trj

  integer :: i, j, k
  integer :: ntmp, wat_idx, mem_comp_idx
  integer :: frame_start, frame_stop, num_components, max_heavy, number_of_frames, number_of_atoms
  integer, allocatable :: component_start(:), num_heavy(:), occupancy(:,:)
  real, allocatable :: criteria(:,:), r(:,:,:)
  integer, allocatable :: atom_map(:,:,:)
  real, dimension(3,3) :: box
  character(len=40) :: fname, iname
  character(len=100) :: line
  integer :: is_water, occ_atom, occ_atom_type, num_mem_components
  real, dimension(3) :: coord = 0.0
  real, allocatable :: closest_hydr(:,:,:)
  character(len=100) :: filename
  integer, allocatable :: atoms_per_component(:)
  character(len=40) :: mapfile





  ! Open Log File
  open(20,file='hydration_shell.log')
  write(20,*) "Beginning hydration calculation"

  
  ! Read Input File
  ! Note - allocates num_heavy, component_start, criteria
  call Alt_Hyd_Input(mapfile, frame_start, frame_stop, fname, iname &
                          , num_components, num_heavy, component_start &
                          , atoms_per_component &
                          , is_water, criteria)

  ! Open Output File
  do i=1, num_components
    write(filename,'(a,i0,a)') 'hydration_shell_',i,'.dat'
    open(100+i,file=trim(filename),status='replace')
  enddo
  
  max_heavy = maxval(num_heavy)
  write(*,*) "The maximum number of heavy atoms per component is ", max_heavy

  ! Call Allocations
  allocate(r(num_components,max_heavy,3))
  allocate(atom_map(num_components,max_heavy,2))
  allocate(closest_hydr(num_components,num_heavy(is_water),2))
  allocate(occupancy(num_components,4))

  ! Zero Arrays
  r = 0.0
  atom_map = 0
  closest_hydr = 0.0
  occupancy = 0

  ! Call Generate_Atom_Map
  call Generate_Atom_Map(mapfile, atoms_per_component, component_start, num_heavy, atom_map)

  ! Open Trajectory
  call trj%open(trim(fname),trim(iname))
  number_of_frames = trj%nframes
  number_of_atoms  = trj%natoms()


  write(20,*) "There are ", number_of_frames," frames."
  write(20,*) "There are ", number_of_atoms," atoms."

  ! Print Frame Information
  do i=1, frame_start-1
    ntmp = trj%read_next(1)
    if (mod(i,100) == 0) write(*,*) "Cycled through frame ", i
  enddo

  !$OMP PARALLEL default(shared)
  frames: do i=frame_start,frame_stop
    !$OMP SINGLE
    write(*,*) "Reached frame ", i
    ! Read Frame and store atoms to ntmp
    ntmp = trj%read_next(1)
    ! Read the box dimensions
    box  = trj%box(1)
    ! Assign membrane Atoms using the atom map
    ! The atom_map is a DD array that points to the atom index and type
    ! Should be 1 entry for each heavy atom (including water O) in order of structure
    comps: do j=1, num_components
      heavs: do k=1, num_heavy(j)
        ! Check that the atom_map is properly initialized
        if (atom_map(j,k,1) == 0) then
          write(*,*) "Error: atom_map is not properly initialized"
          write(*,*) j, k, atom_map(j,k,1)
          stop
        endif
        coord(:) = trj%x(1,atom_map(j,k,1))
        r(j,k,:) = coord(:)
      enddo heavs
    enddo comps
    !$OMP END SINGLE

    !$OMP BARRIER
    
    occupancy = 0
    closest_hydr=0
    !$OMP DO reduction(+:occupancy) private(mem_comp_idx, wat_idx, occ_atom_type)
    waters: do wat_idx=1, num_heavy(is_water)
      ! Loop over lipid and lipid-like atoms to find closest membrane heavy atom
      mem_comps: do mem_comp_idx=1, num_components
        if ( mem_comp_idx /= is_water ) then
          if ( num_heavy(mem_comp_idx) > 0 ) then
            call Find_Hydr(r(is_water,wat_idx,:), r(mem_comp_idx,:,:), num_heavy(mem_comp_idx), mem_comp_idx, box &
                          , atom_map(mem_comp_idx,:,:), criteria, closest_hydr(mem_comp_idx,wat_idx,:))
          endif
          ! Update the occupancy counts
          occ_atom = int(closest_hydr(mem_comp_idx,wat_idx,1))
          if ( occ_atom > 0 ) then
            occ_atom_type = atom_map(mem_comp_idx,occ_atom,2)
            else
            occ_atom_type = 0
          endif
          if (occ_atom_type > 0) then
            occupancy(mem_comp_idx,occ_atom_type) = occupancy(mem_comp_idx,occ_atom_type) + 1
          endif
        endif
      enddo mem_comps

    enddo waters
    !$OMP END DO

    !$OMP SINGLE
    ! Write the occupancy of each type to the log file
    line = ''
    write(*,*) size(occupancy)
    write(20,*) i, ((occupancy(j,k), k=1,4), j=1,num_components)


    ! Write the hydration files
    do mem_comp_idx=1, num_components
      write(*,*) sum(occupancy(mem_comp_idx,:))
      do wat_idx=1, num_heavy(is_water)
        write(100+mem_comp_idx,*) i, closest_hydr(mem_comp_idx,wat_idx,1), closest_hydr(mem_comp_idx,wat_idx,2)
      enddo
    enddo 
    !$OMP END SINGLE
    !$OMP BARRIER

  enddo frames
  !$OMP END PARALLEL

  write(20,*) "Calculation complete"
  close(20)
  close(21)


end program hyd_shell


