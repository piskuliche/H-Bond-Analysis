module hba_module
    implicit none
    public
contains

    subroutine read_hb_input(mapfile, frame_start, frame_stop, fname, iname &
                            , num_components, num_mol, num_acceptors, component_start &
                            , atoms_per_component, is_water, do_water, criteria)

        implicit none

        integer :: i
        integer, intent(out) :: frame_start, frame_stop, num_components, is_water, do_water
        integer, allocatable, intent(out) :: num_mol(:), num_acceptors(:), component_start(:), atoms_per_component(:)
        integer, allocatable :: atoms_per_mol(:), num_acc_sites_per_mol(:), component_label(:)
        integer :: check_if_water, atom_count

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

    subroutine find_h_bonds(r_don, r_acc, num_donors, num_acceptors, criteria, box, atom_map, hydrogen_bonds)
        use myfuncs
        implicit none

        integer :: i, j, hi, oindex
        integer :: num_donors, num_acceptors

        real, intent(in) :: r_don(:,:), r_acc(:,:)
        real, intent(in) :: criteria(:), box(:,:)
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

end module hba_module