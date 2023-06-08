program h_bonding

    use gmxfort_trajectory
    use myfuncs
    use share_routine
    use hba_module

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

    
    integer, allocatable :: atom_map(:,:,:)
    real, dimension(3) :: coord = 0.0
    real, dimension(chunk_size, 3, 3) :: box
    integer, allocatable :: hydrogen_bonds(:,:,:)
    real, allocatable :: r(:,:,:,:)

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
    r = 0.0
    atom_map=0
    hydrogen_bonds = 0

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
    num_frames = frame_stop - frame_start + 1
    nchunks = int(ceiling(real(num_frames) / real(chunk_size)))

    ! Chunking loop
    chunks: do chunk_idx=1, nchunks
        write(*,*) "Reached chunk ", chunk_idx

        ! ****************** CHUNK SETUP ******************
        ! Zero Arrays
        hydrogen_bonds = 0
        box = 0.0

        ! Set the chunk loop start and stop
        chunk_start = (chunk_idx-1)*chunk_size
        chunk_stop = chunk_size
        chunk_stop = min(chunk_size, num_frames - chunk_start)

        ! ****************** READ TRAJECTORY ******************
        write(*,*) "Reading frames ", chunk_start + 1, " to ", chunk_start + chunk_stop
        read_frames_loop: do fr_idx=1, chunk_stop
            ! Read frames
            write(*,*) fr_idx, chunk_start + fr_idx
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
        write(*,*) "Calculating hydrogen bonds"
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(fr_idx, comp_idx)
        frames_loop: do fr_idx=1, chunk_stop
            if ( do_water == 1 ) then
                ! Do water analysis - note this is separated because it is expensive!
                call find_h_bonds(r(fr_idx, is_water, :, :), r(fr_idx, is_water, :,:), num_mol(is_water), num_mol(is_water) &
                                , criteria(is_water,:), box(fr_idx, :, :) &
                                , atom_map(is_water, :, :), hydrogen_bonds(fr_idx, :, :))
            else
                ! Do non-water analysis
                comp_loop: do comp_idx=1, num_components
                    if ( comp_idx /= is_water) then
                        call find_h_bonds(r(fr_idx,is_water,:,:), r(fr_idx, comp_idx, :, :), num_mol(is_water) &
                                        , num_acceptors(comp_idx), criteria(comp_idx,:), box(fr_idx,:,:) &
                                        , atom_map(comp_idx, :,:), hydrogen_bonds(fr_idx, :, :))
                    endif 
                end do comp_loop
            endif
        end do frames_loop
        !$OMP END PARALLEL DO

        ! ****************** OUTPUT DATA ******************
        write(*,*) "Writing data to file"
        write_frames_loop: do fr_idx=1, chunk_stop
            ! DO SOMETHING
            wat: do i=1, num_mol(is_water)
                write(30,*) i, hydrogen_bonds(fr_idx, i, 1), hydrogen_bonds(fr_idx, i, 2)
            end do wat
        end do write_frames_loop

    end do chunks
    close(20)
    close(30)

end program h_bonding

