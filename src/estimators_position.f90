module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,init_path, whichestim, splitbead

  private

  logical, allocatable:: whichestim(:)
  integer:: splitbead

contains

  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof

    nestim=3

    ntime= NMC
    allocate(whichestim(nestim))

    return
  end subroutine init_estimators

  subroutine finalize_estimators()
    implicit none
    deallocate(whichestim)
  end subroutine finalize_estimators

  subroutine estimator(x,p,tcfval)
    implicit none
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision, intent(out):: tcfval(:)
    double precision:: s
    integer:: i,j,k, idof

       tcfval(:)=0.0d0

       do i=1,n
          tcfval(1)=tcfval(1)+x(i,1,1)
          tcfval(2)=tcfval(2)+x(i,1,1)**2
          tcfval(3)=tcfval(3)+x(i,1,1)**3
       end do
       tcfval(:)=tcfval(:)/dble(N)

    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors,weight)
    double precision, intent(out):: x(:,:,:), p(:,:,:),weight,factors(:)
    double precision::              potvals, ringpot, potdiff


    call harmonicsampling(x, p, potvals)

    ringpot=UM(x)
    potdiff= (ringpot- potvals) 
    weight= exp(-betan*potdiff)
    call estimator(x,p,factors)


    return
  end subroutine init_path


end module estimators
