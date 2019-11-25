module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,transfreqs, init_path, whichestim,&
       transition, normalvec, hess, splitbead

  private

  logical, allocatable:: whichestim(:)
  integer:: splitbead
  double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:)
  double precision, allocatable:: lambda(:,:), hess(:,:)

contains

  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof
    !TODO: estimatormod- implement list of estimators

    nestim=3

    ntime= NMC
    allocate(whichestim(nestim))

    allocate(lambda(n,ndof))
    do i=1,n
       do j=1,ndof
          lambda(i,j)= sqrt(transfreqs(j) + lam(i)**2)
       end do
    end do

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
    double precision:: centroidxi, centroidpi, centroidxj, centroidpj
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
    double precision, allocatable:: tempx(:),pk(:),rk(:,:),rp(:), grad(:,:,:), centroidx(:,:)
    double precision::              stdev, potvals, ringpot, potdiff, rpvals
    double precision::              prob, stdevharm, energy, totp, potsum
    integer::                       i,j,k, idof, imin(1)

    allocate(rp(n),tempx(totdof), pk(n))
    allocate(rk(n,ndof))
    stdev= 1.0d0/sqrt(betan)
    stdevharm= 1.0d0/sqrt(beta)
    x(:,:,:)= 0.0d0
    p(:,:,:)= 0.0d0
    potvals=0.0d0
    rpvals=0.0d0

    do i=1,ndim
       do j=1,natom
          idof= calcidof(i,j)
          if (n .gt. 1) then
             !---------------------------
             !generate random numbers for chain
             errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,rp,0.0d0,stdev)!ringpolymer
             !---------------------------
             !scale them to have the right standard deviation
             do k=1, n
                   rp(k)= rp(k)/(lambda(k,idof)*sqrt(mass(j)))
                potvals=potvals + 0.5d0*mass(j)*lambda(k,idof)**2*rp(k)**2
             end do
             !---------------------------
             !transform to site-local cartesian
             if (ring) then
                if (use_fft) then
                   call ring_transform_backward_nr(rp, rk(:,idof))
                else
                   call ring_transform_backward(rp, rk(:,idof))
                end if
             else
                call linear_transform_backward(rp, rk(:,idof), idof)
             end if
          else
             rk(1,idof)=0.0d0
          end if
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,p(:,i,j),0.0d0,stdev)!momenta
          p(:,i,j)= p(:,i,j)*sqrt(mass(j))
       end do
    end do
    !---------------------------
    !transform to global cartesian coordinates
    do k=1,n
       x(k,:,:)= reshape(matmul(transpose(hess),rk(k,:)), (/ndim,natom/))
    end do
    ! errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,ndof,tempx,0.0d0,stdev)
    ! do j=1,ndof
    !    tempx(j)= tempx(j)/sqrt(transfreqs(j))!needs mass scaling, transfreqs are actual freqs
    !    potvals=potvals + 0.5d0*transfreqs(j)*tempx(j)**2 !needs mass scaling, transfreqs are actual freqs
    ! end do
    ! if (ndof.gt.1) then
    !    x(1,:,:)= reshape(matmul(transpose(hess),tempx), (/ndim,natom/))
    !    do i=1,n
    !       do j=1,ndim
    !          do k=1,natom
    !             idof=calcidof(j,k)
    !             x(i,j,k)= x(i,j,k)/sqrt(mass(k)) + rk(i,idof) !does rk not need mass scaling too?
    !          end do
    !       end do
    !    end do
    ! else
    !    do i=1,n
    !       x(i,1,1)= tempx(1)/sqrt(mass(1)) + rk(i,1)
    !    end do
    ! end if
    
    ringpot=UM(x)
    ! do i=1,n
    !    ringpot= ringpot+ pot(x(i,:,:))
    ! end do
    ! allocate(centroidx(ndim,natom))
    ! centroidx(:,:)=0.0d0
    ! do i=1,ndim
    !    do j=1,natom
    !       centroidx(i,j)= centroid(x(:,i,j))
    !    end do
    ! end do
    ! ringpot=pot(centroidx)
    ! deallocate(centroidx)
    potdiff= (ringpot- potvals) 
    weight= exp(-betan*potdiff) !min(exp(-betan*potdiff),1.0d0)
    !work out initial current

    ! write(*,*) potvals, ringpot, potdiff, weight

    call estimator(x,p,factors)

    deallocate(tempx, rp, rk, pk)
    return
  end subroutine init_path


end module estimators
