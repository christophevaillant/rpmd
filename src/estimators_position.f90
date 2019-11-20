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
       do i=1, n-splitbead,1
          s=1.0d0
          do j=0, splitbead-1
             s=s*x(i+j,1,1)
          end do
          tcfval(1)=tcfval(1) + s/dble(N)
          tcfval(2)=tcfval(2) + x(i,1,1)**splitbead/dble(N)
          tcfval(3)= tcfval(3) + x(i,1,1)/dble(N)
       end do
       tcfval(3)= tcfval(3)**splitbead
       tcfval(1)=tcfval(1) + x(n-1,1,1)*x(n,1,1)*x(1,1,1)/dble(N)
       tcfval(1)=tcfval(1) + x(n,1,1)*x(1,1,1)*x(2,1,1)/dble(N)
    ! if (splitbead.gt.1) then
    !    tcfval(1)=0.0d0
    !    do i=1, n-splitbead,splitbead
    !       s=1.0d0
    !       do j=0, splitbead-1
    !          s=s*x(i+splitbead,1,1)
    !       end do
    !       tcfval(1)=tcfval(1) + s/dble(N)
    !    end do
    ! else
    !    do i=1,nestim
    !       tcfval(i)= centroid(x(:,1,1)**3)
    !    end do
    ! end if

    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors,weight)
    double precision, intent(inout):: x(:,:,:), p(:,:,:),weight,factors(:)
    double precision, allocatable:: tempx(:),rk(:,:),rp(:), grad(:,:,:), centroidx(:,:)
    double precision::              stdev, potvals, ringpot, potdiff
    double precision::              prob, stdevharm, energy, totp
    integer::                       i,j,k, idof, imin(1)

    allocate(rp(n),tempx(ndof))
    allocate(rk(n,ndof))
    stdev= 1.0d0/sqrt(betan)
    x(:,:,:)= 0.0d0
    p(:,:,:)= 0.0d0

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
                if (k.eq. 1) then
                   rp(k)= 0.0d0
                else
                   rp(k)= rp(k)/lam(k)
                end if
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
    stdevharm= 1.0d0/sqrt(beta)
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,ndof,tempx,0.0d0,stdevharm)
    tempx(:)= tempx(:)/sqrt(transfreqs(:))!needs mass scaling, transfreqs are actual freqs
    potvals=0.0d0
    do k=1, ndof
       potvals=potvals + 0.5d0*transfreqs(k)*tempx(k)**2 !needs mass scaling, transfreqs are actual freqs
    end do
    if (ndof.gt.1) then
       x(1,:,:)= reshape(matmul(transpose(hess),tempx), (/ndim,natom/))
       do i=1,n
          do j=1,ndim
             do k=1,natom
                idof=calcidof(j,k)
                x(i,j,k)= x(1,j,k)/sqrt(mass(k)) + rk(i,idof) !does rk not need mass scaling too?
             end do
          end do
       end do
    else
       x(1,:,:)= tempx(1)
       do i=1,n
          x(i,1,1)= x(1,1,1)/sqrt(mass(1)) + rk(i,1)
       end do
    end if
    
    ! ringpot=0.0d0
    ! do i=1,n
    !    ringpot= ringpot+ pot(x(i,:,:))
    ! end do
    allocate(centroidx(ndim,natom))
    centroidx(:,:)=0.0d0
    do i=1,ndim
       do j=1,natom
          centroidx(i,j)= centroid(x(:,i,j))
       end do
    end do
    ringpot=pot(centroidx)
    deallocate(centroidx)
    potdiff= (ringpot- potvals)
    weight= exp(-beta*potdiff)
    !work out initial current
    call estimator(x,p,factors)
    deallocate(tempx, rp, rk)
    return
  end subroutine init_path


end module estimators
