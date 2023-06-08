program h_bonding

    use gmxfort_trajectory
    use myfuncs
    use share_routine

    implicit none

    type (Trajectory) :: trj

    ! Loop indexes
    integer :: i, j, k
    integer :: chunk_idx, chunk_start, chunk_stop
    integer :: fr_idx, comp_idx, acc_idx

    ! Input Args
    integer, parameter :: chunk_size=100
    integer :: frame_start, frame_stop
    integer :: num_components, atom_count, is_water, do_water
    integer, allocatable :: num_acceptors(:), component_start(:), atoms_per_component(:), num_mol(:)
    real, allocatable :: criteria(:,:)
    character(len=40) :: mapfile, fname, iname

    ! Main Program
    integer :: number_of_frames, number_of_atoms, ntmp, num_frames
    integer :: nchunks, num_donors, max_acceptors
    integer :: check_if_water
    
    integer, allocatable :: atom_map(:,:,:)
    real, dimension(3) :: coord = 0.0
    real, dimension(chunk_size, 3, 3) :: box
    real, allocatable :: hydrogen_bonds(:,:,:), r(:,:,:,:)

    open(20, file='hydrogen_bonding.log')
    write(20,*) 'Beginning hydrogen bonding analysis'

    ! ****************************************************************************************************
    ! ****************************************************************************************************
    ! Read input file
    call read_hb_input(mapfile, frame_start, frame_stop, fname, iname &
                    , num_components, num_mol, num_acceptors, component_start &
                    , atoms_per_component, is_water, do_water, criteria)

    ! Open Output Files
    open(30, file='all_hydrogen_bonds.dat')
    
    max_acceptors = maxval(num_acceptors)
    num_donors = num_mol(is_water)
    ! Allocate Arrays
    allocate(r(chunk_size, num_components, max_acceptors, 3))
    allocate(atom_map(num_components, max_acceptors, 2))
    allocate(hydrogen_bonds(chunk_size, num_donors, 2))

    ! Zero Arrays

    ! Map Donors and Acceptors
    call Generate_Atom_Map(mapfile, atoms_per_component, component_start, num_acceptors, atom_map)

    ! Open Trajectory
    call trj%open(trim(fname), trim(iname))
    number_of_frames = trj%nframes
    number_of_atoms = trj%natoms()
    ! Print frame information for early frames
    write(20,*) "There are ", number_of_frames, " frames in the trajectory"
    write(20,*) "There are ", number_of_atoms, " atoms in the trajectory"

    ! Skip early frames, and print information
    do i=1, frame_start-1
        ntmp = trj%read_next(1)
        if (mod(i,100) == 0) write(*,*) "Cycled through frame ", i
    end do

    ! ****************************************************************************************************
    ! ****************************************************************************************************

    ! Setup Chunking
    num_frames = frame_stop - frame_start
    nchunks = ceiling(num_frames / chunk_size)

    ! Chunking loop
    chunks: do chunk_idx=1, nchunks
        write(*,*) "Reached chunk ", chunk_idx

        ! ****************** CHUNK SETUP ******************
        ! Zero Arrays
        hydrogen_bonds = 0
        box = 0.0

        ! Set the chunk loop start and stop
        chunk_start = (chunk_idx-1)*chunk_size + 1
        chunk_stop = chunk_idx*chunk_size

        if ( chunk_stop + frame_start > frame_stop ) then
            chunk_stop = frame_stop-chunk_start
        endif

        ! Error Test for problem with chunking
        if ( chunk_stop < chunk_start ) then
            write(*,*) "Error: improper chunking detected"
            write(*,*) "chunk_start: ", chunk_start
            write(*,*) "chunk_stop: ", chunk_stop
            write(*,*) "chunk_idx: ", chunk_idx
            stop 
        endif

        ! ****************** READ TRAJECTORY ******************
        read_frames_loop: do fr_idx=chunk_start, chunk_stop
            ! Read frames
            ntmp = trj%read_next(1)
            box(fr_idx,:,:) = trj%box(1)
            comps_loop: do comp_idx=1, num_components
                accs_loop: do acc_idx=1, num_acceptors(comp_idx)
                    if (atom_map(comp_idx, acc_idx, 1) == 0) then
                        write(*,*) "Error: atom_map is not properly initialized"
                        write(*,*) comp_idx, acc_idx, atom_map(comp_idx, acc_idx, 1)
                        stop
                    endif
                    coord(:) = trj%x(1,atom_map(comp_idx, acc_idx, 1))
                    r(fr_idx,comp_idx, acc_idx, :) = coord(:)
                end do accs_loop
            end do comps_loop
        end do read_frames_loop

        ! ****************** CALCULATION OF HBONDS ******************

        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(fr_idx, comp_idx)
        frames_loop: do fr_idx=chunk_start, chunk_stop
            if ( do_water == 1 ) then
                ! Do water analysis - note this is separated because it is expensive!
                call find_h_bonds(r(fr_idx, is_water, :, :), r(fr_idx, is_water, :,:), num_mol(is_water), num_mol(is_water) &
                                , criteria(is_water,:) &
                                , atom_map(is_water, :, :), hydrogen_bonds(fr_idx, :, :))
            else
                ! Do non-water analysis
                comp_loop: do comp_idx=1, num_components
                    if ( comp_idx /= is_water) then
                        call find_h_bonds(r(fr_idx,is_water,:,:), r(fr_idx, comp_idx, :, :), num_mol(is_water) &
                                        , num_acceptors(comp_idx), criteria(comp_idx,:) &
                                        , atom_map(comp_idx, :,:), hydrogen_bonds(fr_idx, :, :))
                    endif 
                end do comp_loop
            endif
        end do frames_loop
        !$OMP END PARALLEL DO

        ! ****************** OUTPUT DATA ******************
        write_frames_loop: do fr_idx=chunk_start, chunk_stop
            ! DO SOMETHING
            wat: do i=1, num_mol(is_water)
                write(30,*) i, hydrogen_bonds(fr_idx, i, 1), hydrogen_bonds(fr_idx, i, 2)
            end do wat
        end do write_frames_loop

    end do chunks
    close(20)
    close(30)

end program h_bonding

subroutine read_hb_input(mapfile, frame_start, frame_stop, fname, iname &
                        , num_components, num_mol, num_acceptors, component_start &
                        , atoms_per_component, is_water, do_water, criteria)

    implicit none

    integer :: i
    integer, intent(out) :: frame_start, frame_stop, num_components, atom_count, is_water, do_water
    integer, allocatable, intent(out) :: num_mol(:), num_acceptors(:), component_start(:), atoms_per_component(:)
    integer, allocatable :: atoms_per_mol(:), num_acc_sites_per_mol(:), component_label(:)

    real, allocatable :: criteria(:,:)

    character(len=40) :: mapfile, fname, iname



    open(10, file='hbonding.in', status='old')
    read(10,*)
    read(10,*) fname, iname, mapfile
    read(10,*)
    read(10,*) frame_start, frame_stop
    read(10,*)
    read(10,*) num_components, do_water
    read(10,*)

    allocate(num_mol(num_components))
    allocate(atoms_per_mol(num_components))
    allocate(num_acc_sites_per_mol(num_components))
    allocate(num_acceptors(num_components))
    allocate(atoms_per_component(num_components))
    allocate(component_start(num_components))
    allocate(criteria(num_components,3))

    num_mol = 0; atoms_per_mol = 0; num_acc_sites_per_mol = 0
    num_acceptors = 0; atoms_per_component = 0; component_start = 0
    criteria = 0.0
  

    atom_count = 0
    is_water = -1
    do i=1, num_components
        check_if_water = -1
        read(10,*) component_label(i), num_mol(i), atoms_per_mol(i), num_acc_sites_per_mol(i)
        read(10,*) is_water

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
        num_acceptors(i) = num_acc_sites_per_mol(i)*num_mol(i)
        ! Number of atoms per component is the number of molecules times the number of atoms per molecule
        atoms_per_component(i) = num_mol(i)*atoms_per_mol(i)
        ! Total number of atoms is the sum of the atoms per component
        atom_count = atom_count + atoms_per_component(i)
    end do
    read(10,*)

    ! Read in the h-bond criteria
    do i=1, num_components
        ! OX distance, HX distance, OX-HX angle
        read(10,*) criteria(i,1), criteria(i,2), criteria(i,3)
        criteria(i,1) = criteria(i,1)**2
        criteria(i,2) = criteria(i,2)**2
    end do

    close(10)

    write(*,*) "Finished reading input file"

end subroutine read_hb_input

subroutine find_h_bonds(r_don, r_acc, criteria, num_donors, num_acceptors, box, atom_map, hydrogen_bonds)
    use myfuncs
    implicit none

    integer :: i, j, hi, oindex
    integer :: num_donors, num_acceptors

    real, intent(in) :: r_don(:,:), r_acc(:,:)
    real, intent(in) :: criteria(:,:), box(:,:)
    real :: roxsq, rhxsq, angle
    real, dimension(num_donors, 2), intent(out) :: hydrogen_bonds
    real, intent(in) :: atom_map(:,:)

    ! Loop over oxygens
    donors: do i=1, num_donors
        oindex = (i-1)*3 + 1
        acceptors: do j=1, num_acceptors
            roxsq = distance2(r_don(oindex,:), r_acc(i,:), box)
            if ( roxsq < criteria(1) ) then
                ! Loop over hydrogens
                hatoms: do hi=1,2
                    rhxsq = distance2(r_don(oindex+hi,:), r_acc(i,:), box)
                    if ( rhxsq < criteria(2) ) then
                        ! Calculate angle
                        angle = bond_angle(r_don(oindex+hi,:), r_don(oindex,:), r_acc(j,:), box)
                        if ( angle < criteria(3) ) then
                            ! COUNT DATA
                            hydrogen_bonds(i, hi) = atom_map(j,2)
                            ! We can exit, because a single water can only donate one hydrogen bond to a single acceptor
                            exit hatoms
                        endif
                    endif
                enddo hatoms
            endif
        end do acceptors
        ! Loop over hydrogens
    end do donors

end subroutine find_h_bonds