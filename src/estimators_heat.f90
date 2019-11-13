module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,transfreqs, init_path, whichestim,&
       transition, normalvec, hess, normalization, convection

  private

  logical, allocatable:: whichestim(:)
  logical:: convection
  double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:)
  double precision, allocatable:: lambda(:,:), hess(:,:)

contains

  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof
    !TODO: estimatormod- implement list of estimators

    if (outputfbar) open(200,file="fbar.dat")

    ntime= NMC/Noutput
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
    double precision:: s, energy, prob, stdev
    integer:: i,j,k, idof

    tcfval(1)=0.0d0
    do i=1,natom
       energy=0.0d0
       do j=1,n
          ! if (i .lt. natom) energy= energy + &
          !      (x(j,1,i+1) - x(j,1,i))*interforce(x,i,j)*p(j,1,i)/mass(i)
          ! if (i .gt. 1) energy= energy + &
          !      (x(j,1,i) - x(j,1,i-1))*interforce(x,i-1,j)*p(j,1,i)/mass(i)
          if (i .lt. natom) energy= energy + &
               0.5d0*(x(j,1,i+1) - x(j,1,i))*interforce(x,i,j)*(p(j,1,i) + p(j,1,i+1))/mass(i)
          if (i .gt. 1) energy= energy + &
               0.5d0*(x(j,1,i) - x(j,1,i-1))*interforce(x,i-1,j)*(p(j,1,i-1) + p(j,1,i))/mass(i)
          if (convection) energy= energy+ p(j,1,i)*siteenergy(p,x,i,j)/mass(i)
       end do
       tcfval(1)= tcfval(1)+energy/dble(N)
    end do

    tcfval(2)= tcfval(1)

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
    ! do k=1,n
    !    do i=1,ndim
    !       totp=0.0d0
    !       do j=1,natom
    !          totp=totp+ p(k,i,j)
    !       end do
    !       p(k,i,:)= p(k,i,:) - totp
    !    end do
    ! end do
    !---------------------------
    !transform to global cartesian coordinates
    stdevharm= 1.0d0/sqrt(beta)
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,ndof,tempx,0.0d0,stdevharm)
    tempx(:)= tempx(:)/sqrt(transfreqs(:))
    potvals=0.0d0
    do k=1, ndof
       potvals=potvals + 0.5d0*transfreqs(k)*tempx(k)**2 
    end do
    x(1,:,:)= reshape(matmul(transpose(hess),tempx), (/ndim,natom/))
    do i=1,n
       x(i,:,:)= x(1,:,:)
    end do
    do i=1,ndim
       do j=1,natom
          idof= calcidof(i,j)
          x(:,i,j)= x(:,i,j)/sqrt(mass(j)) + lattice(i,j) + rk(:,idof)
       end do
    end do
    
    ! ringpot= 0.0d0
    ! do k=1,n
    !    ringpot= ringpot+ pot(x(k,:,:))
    ! end do
    allocate(centroidx(ndim,natom))
    do i=1,ndim
       do j=1,natom
          centroidx(i,j)= centroid(x(:,i,j))
       end do
    end do
    ringpot= pot(centroidx)
    potdiff= (ringpot- potvals)
    weight= exp(-beta*potdiff)
    !work out initial current
    call estimator(x,p,factors)
    factors(2)=1.0d0
    
    deallocate(tempx, rp, rk, centroidx)
    return
  end subroutine init_path

  function normalization()
    implicit none
    double precision:: normalization
    integer:: i
    
    normalization=1.0d0 !exp(-beta*V0)
    do i=1, ndof
       normalization=normalization*harmonicpart(transfreqs(i)*beta**2)
    end do

    return
  end function normalization

  function harmonicpart(omegasq)
    implicit none
    double precision:: harmonicpart
    double precision, intent(in):: omegasq
    integer:: i

    harmonicpart= 1.0d0
    do i=1,n
       harmonicpart= harmonicpart/sqrt(4.0d0*sin(dble(i)*pi/dble(n))**2 + omegasq/(n**2))
    end do

    return
  end function harmonicpart

end module estimators
