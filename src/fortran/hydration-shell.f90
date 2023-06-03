program hyd_shell
  use hydr_module
  use gmxfort_trajectory
  use myfuncs
  implicit none

  type (Trajectory) :: trj

  integer :: i, j, k
  integer :: ntmp, wat_idx, mem_comp_idx
  integer :: frame_start, frame_stop, num_components, max_heavy, number_of_frames, number_of_atoms
  integer, allocatable :: component_start(:), num_heavy(:), occupancy(:)
  real, allocatable :: criteria(:), r(:,:,:)
  integer, allocatable :: atom_map(:,:,:)
  real, dimension(3,3) :: box
  character(len=40) :: fname, iname
  character(len=100) :: line
  integer :: is_water, occ_atom_type, num_mem_components
  real, dimension(3) :: coord = 0.0
  real, allocatable :: closest_hydr(:,:,:)
  character(len=100) :: filename
  integer, allocatable :: atoms_per_component(:)





  ! Open Log File
  open(20,file='hydration_shell.log')
  write(20,*) "Beginning hydration calculation"

  
  ! Read Input File
  ! Note - allocates num_heavy, component_start, criteria
  call Alt_Hyd_Input(frame_start, frame_stop, fname, iname &
                          , num_components, num_heavy, component_start &
                          , atoms_per_component &
                          , is_water, criteria)

  ! Open Output File
  do i=1, num_components
    write(filename,'(a,i0,a)') 'hydration_shell_',i,'.dat'
    open(100+i,file=trim(filename),status='replace')
  enddo
  
  max_heavy = maxval(num_heavy)

  ! Call Allocations
  allocate(r(num_components,max_heavy,3))
  allocate(atom_map(num_components,max_heavy,2))
  allocate(closest_hydr(num_components,num_heavy(is_water),2))
  allocate(occupancy(size(criteria)))
  ! Zero Arrays
  r = 0.0
  atom_map = 0
  closest_hydr = 0.0
  occupancy = 0


  ! Maps the heavy atoms to for the trajectory
  ! TODO: FIX THIS
  call Generate_Atom_Map(atoms_per_component, component_start, num_heavy, atom_map)

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

  !$OMP PARALLEL
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
    comps: do j=1,num_components
      heavs: do k=1,num_heavy(j)
        coord(:) = trj%x(1,atom_map(j,k,1))
        r(j,k,:) = coord(:)
      enddo heavs
    enddo comps
    !$OMP END SINGLE

    !$OMP BARRIER
    
    occupancy = 0
    !$OMP DO PRIVATE(closest_hydr,mem_comp_idx) reduction(+:occupancy) schedule(DYNAMIC)
    waters: do wat_idx=1, num_heavy(is_water)
      ! Loop over lipid and lipid-like atoms to find closest membrane heavy atom
      mem_comps: do mem_comp_idx=1, num_mem_components
        if ( mem_comp_idx /= is_water ) then
          closest_hydr=0
          call Find_Hydr(r(is_water,wat_idx,:), r(mem_comp_idx,:,:), num_heavy(mem_comp_idx), box &
                          , atom_map(mem_comp_idx,:,:), criteria, closest_hydr(mem_comp_idx,wat_idx,:))
          ! Update the occupancy counts
          occ_atom_type = closest_hydr(mem_comp_idx,wat_idx,2)
          occupancy(occ_atom_type) = occupancy(occ_atom_type) + 1
        endif
      enddo mem_comps

    enddo waters
    !$OMP END DO

    !$OMP SINGLE
    ! Write the occupancy of each type to the log file
    line = ''
    do j=1, size(occupancy)
      write(line,'(a,i0,a,i0)') occupancy(i), ' '
    enddo
    write(20,*) i, line

    ! Write the hydration files
    do mem_comp_idx=1, num_mem_components
      do wat_idx=1, num_heavy(is_water)
        write(100+mem_comp_idx,*) i, closest_hydr(mem_comp_idx,wat_idx,1), closest_hydr(mem_comp_idx,wat_idx,2)
      enddo
    enddo 
    !$OMP END SINGLE

  enddo frames
  !$OMP END PARALLEL

  write(20,*) "Calculation complete"
  close(20)
  close(21)


end program hyd_shell


