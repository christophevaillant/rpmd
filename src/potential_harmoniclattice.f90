module potential
  use general
  implicit none
  double precision::               harm, interharm, V0, gammaheat,spacing
  double precision, allocatable::  lattice(:,:)
  double precision:: betaleft, betaright

  public

contains
  subroutine V_init(iproc)
    integer, intent(in):: iproc
    integer:: i
    namelist /POTDATA/ spacing, harm, interharm, gammaheat

    spacing=1.0d0
    harm=1.0d-2
    interharm=1.0d-2
    gammaheat=1.0d0

    if (iproc.eq.0) read(5, nml=POTDATA)

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

    siteenergy=0.5d0*p(j,1,i)**2/mass(i)
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
    !langevin thermostat step
  subroutine pot_thermostat(x,p)
    implicit none
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision,allocatable:: v(:,:), newp(:), grad(:,:,:)
    double precision:: a1,a2
    integer:: i,j,k, idof

    allocate(newp(ndof))
    a1=exp(-gammaheat*dt)!
    a2= sqrt(1.0d0- a1**2)
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,ndof,newp,0.0d0,1.0d0)
    !middle atoms
    do i=1, ndim
       do j=2, natom-1
          idof= calcidof(i,j)
          do k=1,n
             p(k,i,j)= (a1**2)*p(k,i,j) + &
                  sqrt(mass(j)/beta)*a2*sqrt(1.0+a1**2)*newp(idof)
          end do
       end do
    end do
    do i=1, ndim
       !left reservoir
       j=1
       idof= calcidof(i,j)
       do k=1,n
          p(k,i,j)= (a1**2)*p(k,i,j) + &
               sqrt(mass(j)/betaleft)*a2*sqrt(1.0+a1**2)*newp(idof)
       end do
       !right reservoir
       j=natom
       idof= calcidof(i,j)
       do k=1,n
          p(k,i,j)= (a1**2)*p(k,i,j) + &
               sqrt(mass(j)/betaright)*a2*sqrt(1.0+a1**2)*newp(idof)
       end do
    end do


    deallocate(newp)
  end subroutine pot_thermostat


  function calcpartition()
    implicit none
    double precision:: calcpartition

    calcpartition= sqrt(mass(1)/(2.0d0*pi*beta))
  end function calcpartition

end module potential
