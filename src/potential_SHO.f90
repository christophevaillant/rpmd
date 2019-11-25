module potential
  use general
  implicit none
  double precision::               omega, V0

  public

contains
  subroutine V_init(iproc)
    integer, intent(in):: iproc
    namelist /POTDATA/ omega

    omega=1.0

    if (iproc.eq.0) read(5, nml=POTDATA)

    V0=0.0d0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function pot(x)
    implicit none
    double precision::     pot, x(:,:), q1, q2
    integer::              i,j

    ! pot= 0.5d0*x(1,1)**2 + 0.1d0*x(1,1)**3 + 1.0d-2*x(1,1)**4
    ! pot= 0.25d0*x(1,1)**4
    pot= 0.5d0*mass(1)*omega**2*x(1,1)**2
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

    ! grad(1,1)=x(1,1)**3
    ! grad(1,1)=x(1,1) + 0.3d0*x(1,1)**2 + 4d-2*x(1,1)**3
    grad(1,1)=mass(1)*omega**2*x(1,1)
    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision, intent(out)::   hess(:,:,:,:)
    double precision, intent(in)::  x(:,:)
    double precision:: q1, q2, q3
    integer::              i, j

    ! hess(1,1,1,1)=3.0d0*x(1,1)**2! + omega**2
    ! hess(1,1,1,1)=1.0d0 + 0.6d0*x(1,1) + 1.2d-1*x(1,1)**2
    hess(1,1,1,1)=mass(1)*omega**2

    return
  end subroutine Vdoubleprime


  function calcpartition()
    implicit none
    double precision:: calcpartition

    calcpartition= sqrt(mass(1)/(2.0d0*pi*beta))
  end function calcpartition

  !---------------------------------------------------------------------
    !langevin thermostat step
  subroutine pot_thermostat(x,p)
    implicit none
    double precision, intent(inout):: x(:,:,:), p(:,:,:)

    return
  end subroutine pot_thermostat

end module potential
