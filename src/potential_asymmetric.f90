module potential
  use general
  implicit none
  double precision::               V0, Aeck, Beck, x0

contains
  subroutine V_init()
    Aeck=-18.0d0/pi
    Beck=13.5d0/pi
    x0=8.0d0/sqrt(3.0d0*pi)
    V0=0.5d0*Aeck+ Beck
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function pot(x)
    implicit none
    double precision::     pot, x(:,:)
    integer::              i,j

    pot= Aeck/(1.0d0+ exp(-2.0d0*x(1,1)/x0))
    pot= pot +Beck/cosh(x(1,1)/x0)**2
    
    return
  end function POT

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:)

    grad(:,:)= Aeck/(x0*(1.0d0+ cosh(2.0d0*x(:,:)/x0)))
    grad(:,:)=grad(:,:) -2.0d0*Beck*(tanh(x(:,:)/x0)/cosh(x(:,:)/x0)**2)/x0

    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j

    hess(1,1,1,1)= -Aeck*tanh(x(1,1)/x0)/(x0*cosh(x(1,1)/x0))**2

    hess(1,1,1,1)= hess(1,1,1,1)+ &
         Beck*(4.0d0*(tanh(x(1,1)/x0)/cosh(x(1,1)/x0))**2 -&
         2.0d0/cosh(x(1,1)/x0)**4)/x0**2

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
