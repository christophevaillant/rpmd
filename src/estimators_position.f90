module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,init_path, whichestim, splitbead

  private

  double precision:: H0
  logical, allocatable:: whichestim(:)
  integer:: splitbead

contains

  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof

    if (nonlinear.eq.0) then
       nestim=3
    else
       nestim=6
    end if
    
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
    double precision, allocatable:: lagrangian(:)
    double precision:: lagsum1, lagsum2
    integer:: i,j,k, idof, jindex,kindex

    tcfval(:)=0.0d0

    if (nonlinear.eq.1) then
       do i=1,n
          tcfval(1)=tcfval(1)+x(i,1,1)/dble(N)
          tcfval(2)=tcfval(2)+x(i,1,1)**2/dble(N**2)
          tcfval(3)=tcfval(3)+x(i,1,1)**3/dble(N**3)
          do j=i+1,n
             tcfval(2)=tcfval(2)+2.0d0*x(i,1,1)*x(j,1,1)/dble(N**2)
             tcfval(3)=tcfval(3)+2.0d0*x(i,1,1)*x(j,1,1)**2/dble(N**3)
             do k=j+1,n
                tcfval(3)=tcfval(3)+6.0d0*x(i,1,1)*x(j,1,1)*x(k,1,1)/dble(N**3)
             end do
          end do
          tcfval(5)=tcfval(5)+x(i,1,1)**2/dble(N)
          tcfval(6)=tcfval(6)+x(i,1,1)**3/dble(N)
       end do
       tcfval(4)=tcfval(1)
    else if (nonlinear.eq.2) then
       do i=1,n
          tcfval(1)=tcfval(1)+x(i,1,1)/dble(N)
          tcfval(2)=tcfval(2)+x(i,1,1)**2/dble(N)
          tcfval(3)=tcfval(3)+x(i,1,1)**3/dble(N)
          do j=i+1,n
             tcfval(2)=tcfval(2)+6.0d0*x(i,1,1)*x(j,1,1)/dble(N**2)
             tcfval(3)=tcfval(3)+6.0d0*x(i,1,1)*x(j,1,1)**2/dble(N**2)
             do k=j+1,n
                tcfval(3)=tcfval(3)+48.0d0*x(i,1,1)*x(j,1,1)*x(k,1,1)/dble(N**3)
             end do
          end do
          tcfval(5)=tcfval(5)+x(i,1,1)**2/dble(N)
          tcfval(6)=tcfval(6)+x(i,1,1)**3/dble(N)
       end do
       tcfval(4)=tcfval(1)
    else
       do i=1,n
          tcfval(1)=tcfval(1)+x(i,1,1)/dble(N)
          tcfval(2)=tcfval(2)+x(i,1,1)**2/dble(N)
          tcfval(3)=tcfval(3)+x(i,1,1)**3/dble(N)
       end do
    end if
    
    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors,weight)
    double precision, intent(out):: x(:,:,:), p(:,:,:),weight,factors(:)
    double precision::              potvals, ringpot, potdiff,throwaway(nestim)
    integer:: i

    call harmonicsampling(x, p, potvals)

    ringpot=UM(x)
    potdiff= (ringpot- potvals) 
    weight= exp(-betan*potdiff)
    ! if (calcphase) then
    !    calcphase=.false.
    !    call estimator(x,p,factors)
    !    calcphase=.true.
    ! else
       call estimator(x,p,factors)
    ! end if
    return
  end subroutine init_path


end module estimators
