module potential
  use general
  implicit none
  double precision::               V0, Vheight, x0

contains
  subroutine V_init()
    Vheight=1.56185d-2
    x0=0.734d0
    V0=Vheight
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function pot(x)
    implicit none
    double precision::     pot, x(:,:)
    integer::              i,j

    !sum only there to make arrays fit
    pot= Vheight/cosh(x(1,1)/x0)**2 !- V0

    return
  end function POT

  !---------------------------------------------------------------------
  subroutine Vprime(x, grad)
    implicit none
    integer::              i,j
    double precision::     grad(:,:), x(:,:)

    grad(:,:)= -2.0d0*Vheight*(tanh(x(:,:)/x0)/cosh(x(:,:)/x0)**2)/x0

    return
  end subroutine Vprime

  !---------------------------------------------------------------------
  subroutine  Vdoubleprime(x,hess)
    implicit none
    double precision::     hess(:,:,:,:), x(:,:), dummy1, eps
    integer::              i, j

    hess(1,1,1,1)= Vheight*(4.0d0*(tanh(x(1,1)/x0)/cosh(x(1,1)/x0))**2 -&
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
