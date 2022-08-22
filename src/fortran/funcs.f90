module myfuncs
  implicit none
  public
contains

    function pbc(a, box)

        implicit none
        real, intent(in) :: a(3), box(3,3)
        real :: pbc(3)
        integer :: I

        pbc = a

        do I = 3, 1, -1
            pbc(1:I) = pbc(1:I) - box(1:I,I) * nint(pbc(I) / box(I,I))
        end do

    end function pbc

    function distance2(a, b, box)

        implicit none
        real :: distance2
        real, intent(in), dimension(3) :: a, b
        real :: c(3)
        real, intent(in) :: box(3,3)

        c = pbc(a - b, box)
        distance2 = dot_product(c, c)

    end function distance2

    function bond_vector(a, b, box)

        implicit none
        real :: bond_vector(3)
        real, intent(in), dimension(3) :: a, b
        real, intent(in) :: box(3,3)

        bond_vector = pbc(a-b, box)

    end function bond_vector


    function magnitude(a)

        implicit none
        real :: magnitude
        real, intent(in) :: a(3)

        magnitude = sqrt(dot_product(a, a))

    end function magnitude


    function bond_angle(a, b, c, box)

        implicit none
        real :: bond_angle
        real, intent(in), dimension(3) :: a, b, c
        real, intent(in) :: box(3,3)
        real, dimension(3) :: bond1, bond2

        bond1 = bond_vector(b, a, box)
        bond2 = bond_vector(b, c, box)

        bond_angle = acos(dot_product(bond1, bond2)/(magnitude(bond1)*magnitude(bond2)))

    end function bond_angle

end module myfuncs
