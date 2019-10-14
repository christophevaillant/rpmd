module potential
  use general
  implicit none
  double precision::               harm, interharm, V0, gammaheat,spacing
  double precision, allocatable::  lattice(:,:)

  public

contains
  subroutine V_init()
    integer:: i
    namelist /POTDATA/ spacing, harm, interharm, gammaheat

    spacing=1.0d0
    harm=1.0d-2
    interharm=1.0d-2
    gammaheat=1.0d0

    read(5, nml=POTDATA)

    allocate(lattice(1,natom))
    do i=1, natom
       lattice(1,i)= i*spacing
    end do
    V0=0.0d0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function pot(x)
    implicit none
    double precision::     pot, x(:,:), q1, q2
    integer::              i,j

    pot=0.0d0
    !lattice potential
    do i=1, natom
       pot= pot + 0.5d0*mass(i)*harm**2*(lattice(1,i)-x(1,i))**2
    end do
    !interatomic potential
    do i=1, natom-1
       q1=x(1,i) - lattice(1,i)
       q2=x(1,i+1) - lattice(1,i+1)
       pot= pot + 0.5d0*mass(i)*interharm**2*(q1-q2)**2
    end do
    pot=pot-V0
    return
  end function POT

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision, intent(in)::   x(:,:)
    double precision, intent(out)::  grad(:,:)
    double precision:: q1, q2, q3

    grad(:,:)=0.0d0
    do i=1, natom
       q1=x(1,i) - lattice(1,i)
       grad(1,i)= mass(i)*harm**2*q1
       if (i.eq.1) then
          q2=x(1,i+1) - lattice(1,i+1)
          grad(1,i)= grad(1,i)+ mass(i)*interharm**2*(q1- q2)
       else if (i.eq.natom) then
          q2=x(1,i-1) - lattice(1,i-1)
          grad(1,i)= grad(1,i)+ mass(i)*interharm**2*(q1-q2)
       else
          q2=x(1,i+1) - lattice(1,i+1)
          q3=x(1,i-1) - lattice(1,i-1)
          grad(1,i)= grad(1,i)+ mass(i)*interharm**2*(2.0d0*q1 - q2 - q3)
       end if
    end do

    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision, intent(out)::   hess(:,:,:,:)
    double precision, intent(in)::  x(:,:)
    double precision:: q1, q2, q3
    integer::              i, j

    hess(:,:,:,:)=0.0d0
    do i=1, natom
       do j=1, natom
          if (i.eq.j) then
             if (i.eq.1 .or. i.eq.natom) then
                hess(1,i,1,j)= mass(i)*harm**2 + mass(i)*interharm**2
             else
                hess(1,i,1,j)= mass(i)*harm**2 + 2.0d0*mass(i)*interharm**2
             end if
          else if (i .eq. j+1 .or. i .eq. j-1) then
             hess(1,i,1,j)= -mass(i)*interharm**2
          end if
       end do
    end do

    return
  end subroutine Vdoubleprime

  !---------------------------------------------------------------------
  !inter-bead force, nearest neighbours, pbc
  function beadforce(x,i,j)
    implicit none
    double precision, intent(in):: x(:,:,:)
    integer, intent(in):: i,j
    double precision:: beadforce, q1, q2

    q1= x(j,1,i)
    if (j .eq. n) then
       q2= x(1,1,i)
    else
       q2= x(j+1,1,i)
    end if
    
    beadforce= mass(i)*(q1-q2)/betan**2

    return
  end function beadforce

  !---------------------------------------------------------------------
  !interatomic force, no pbc
  function interforce(x,i,j)
    implicit none
    double precision, intent(in):: x(:,:,:)
    integer, intent(in):: i,j
    double precision:: interforce, q1, q2

    q1= x(j,1,i) - lattice(1,i)
    if (i .lt. natom) then
       q2= x(j,1,i+1) - lattice(1,i+1)
       interforce= mass(i)*interharm**2*(q1 - q2)
    else
       interforce=0.0d0
    end if
    
    return
  end function interforce

  !---------------------------------------------------------------------
  !site energy
  function siteenergy(p,x,i,j)
    implicit none
    double precision, intent(in):: x(:,:,:),p(:,:,:)
    integer, intent(in):: i,j
    double precision:: siteenergy, q1, q2

    siteenergy=0.0d0
    siteenergy=siteenergy+ 0.5d0*p(j,1,i)**2/mass(i)
    q1= x(j,1,i) - lattice(1,i)
    if (i .lt. natom) then
       q2= x(j,1,i+1) - lattice(1,i+1)
       siteenergy= siteenergy+0.25d0*mass(i)*interharm**2*(q1 - q2)**2
    end if
    if (i .gt. 1) then
       q2= x(j,1,i-1) - lattice(1,i-1)
       siteenergy= siteenergy+0.25d0*mass(i)*interharm**2*(q2 - q1)**2
    end if
    siteenergy=siteenergy+0.5d0*mass(i)*harm**2*q1**2

    return
  end function siteenergy


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
