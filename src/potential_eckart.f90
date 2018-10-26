module potential
  use general
  implicit none
  double precision::               V0, Vheight, x0

contains
  subroutine V_init()
    Vheight=1.0d0
    x0=1.0d0
    return
  end subroutine V_init
  !---------------------------------------------------------------------
  function pot(x)
    implicit none
    double precision::     pot, x(:,:)
    integer::              i,j

    !sum only there to make arrays fit
    pot= Vheight/cosh(x(1,1)/x0)**2

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

end module potential
