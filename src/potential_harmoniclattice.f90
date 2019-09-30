module potential
  use general
  implicit none
  double precision::               harm, interharm
  double precision, allocatable::  lattice(:,:)

  public

contains
  subroutine V_init()
    double precision:: spacing
    integer:: i
    namelist /POTDATA/ spacing, harm, interharm

    spacing=1.0d0
    harm=1.0d-2
    interharm=1.0d-2

    read(5, nml=POTDATA)

    allocate(lattice(1,natom))
    do i=1, natom
       lattice(1,i)= i*spacing
       write(*,*)i, lattice(1,i)
    end do

    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function pot(x)
    implicit none
    double precision::     pot, x(:,:)
    integer::              i,j

    pot=0.0d0
    !lattice potential
    do i=1, natom
       pot= pot + 0.5d0*mass(i)*harm**2*(lattice(1,i)-x(1,i))**2
    end do
    !interatomic potential
    do i=1, natom-1
       pot= pot + 0.5d0*mass(i)*interharm**2*(x(1,i+1)-x(1,i))**2
    end do

    return
  end function POT

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:)

    do i=1, natom
       grad(1,i)= mass(i)*harm**2*(x(1,i)-lattice(1,i))
       if (i.eq.1) then
          grad(1,i)= grad(1,i)+ mass(i)*interharm**2*(x(1,i)- x(1,i+1))
       else if (i.eq.natom) then
          grad(1,i)= grad(1,i)+ mass(i)*interharm**2*(x(1,i)- x(1,i-1))
       else
          grad(1,i)= grad(1,i)+ mass(i)*interharm**2*(2.0d0*x(1,i)- x(1,i-1)- x(1,i+1))
       end if
    end do

    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j

    do i=1, natom
       do j=1, natom
          if (i.eq.j) then
             hess(1,i,1,j)= mass(i)*harm**2 + 2.0d0*mass(i)*interharm**2
          else if (i .eq. j+1 .or. i .eq. j-1) then
             hess(1,i,1,j)= -mass(i)*interharm**2
          end if
       end do
    end do

    return
  end subroutine Vdoubleprime

  !---------------------------------------------------------------------
  subroutine  massweightedhess(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j

    call Vdoubleprime(x,hess)
    hess(1,1,1,1)=hess(1,1,1,1)/mass(1)
    return
  end subroutine Massweightedhess


  function calcpartition()
    implicit none
    double precision:: calcpartition

    calcpartition= sqrt(mass(1)/(2.0d0*pi*beta))
  end function calcpartition

end module potential
