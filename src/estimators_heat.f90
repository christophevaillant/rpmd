module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,transfreqs, init_path, whichestim,&
       betaright, betaleft, width, transition, normalvec, hess, normalization

  private

  logical, allocatable:: whichestim(:)
  double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:)
  double precision, allocatable:: lambda(:,:), hess(:,:)
  double precision:: betaleft, betaright, width

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
    double precision,allocatable:: v(:,:), newp(:), grad(:,:,:)
    double precision:: s, energy, prob, stdev, a1,a2
    integer:: i,j,k, idof

    tcfval(1)=0.0d0
    do i=1,natom
       energy=0.0d0
       do j=1,n
          if (i .lt. natom) energy= energy + &
               0.5d0*(x(j,1,i+1) - x(j,1,i))*interforce(x,i,j)*(p(j,1,i) + p(j,1,i+1))/mass(i)
          energy= energy+ p(j,1,i)*siteenergy(p,x,i,j)/mass(i)
       end do
       tcfval(1)= tcfval(1)+energy/dble(N)
    end do

    !langevin thermostat step
    allocate(newp(ndof))
    a1=exp(-gammaheat*dt)!
    a2= sqrt(1.0d0- a1**2)
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,ndof,newp,0.0d0,1.0d0)
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
       j=1
       idof= calcidof(i,j)
       do k=1,n
          p(k,i,j)= (a1**2)*p(k,i,j) + &
               sqrt(mass(j)/betaleft)*a2*sqrt(1.0+a1**2)*newp(idof)
       end do
       j=natom
       idof= calcidof(i,j)
       do k=1,n
          p(k,i,j)= (a1**2)*p(k,i,j) + &
               sqrt(mass(j)/betaright)*a2*sqrt(1.0+a1**2)*newp(idof)
       end do
    end do


    tcfval(2)= tcfval(1)
    deallocate(newp)

    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors,weight)
    double precision, intent(out):: x(:,:,:), p(:,:,:),weight,factors(:)
    double precision, allocatable:: tempx(:),rk(:,:),rp(:), grad(:,:,:), centroidx(:,:)
    double precision::              stdev, potvals, ringpot, potdiff
    double precision::              prob, stdevharm, energy
    integer::                       i,j,k, idof, imin(1), inn

    allocate(rp(n),tempx(1))
    allocate(rk(n,ndof))
    potvals=0.0d0
    stdev= 1.0d0/sqrt(betan)
    x(:,:,:)= 0.0d0
    p(:,:,:)= 0.0d0

    do idof=1, ndof
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
    end do
    !---------------------------
    !transform to global cartesian coordinates    
    potvals=0.0d0
    do i=1,ndim
       do j=1,natom
          idof= calcidof(i,j)
          stdevharm= 1.0d0/sqrt(beta*transfreqs(idof)*mass(j))
          stdev= 1.0d0/sqrt(betan)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,1,tempx,lattice(i,j),stdevharm)!sites
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,p(:,i,j),0.0d0,stdev)!momenta
          x(:,i,j)= rk(:,idof) + tempx(1)
          do k=1, n
             potvals=potvals + 0.5d0*mass(j)*transfreqs(idof)*(x(k,i,j)-lattice(i,j))**2 
          end do
          p(:,i,j)= p(:,i,j)*sqrt(mass(j))
       end do
    end do

    allocate(centroidx(ndim,natom))
    do i=1,ndim
       do j=1,natom
          centroidx(i,j)= centroid(x(:,i,j))
       end do
    end do
    ringpot= 0.0d0
    do k=1,n
       ringpot= ringpot+ pot(x(k,:,:))
    end do
    potdiff= (ringpot- potvals)
    weight= exp(-betan*potdiff) !min(exp(-betan*potdiff), 1.0d0) !

    !work out initial current
    factors(1)=0.0d0
    do i=1,natom
       energy=0.0d0
       do j=1,n
          if (i .lt. natom) energy= energy + &
               0.5d0*(x(j,1,i+1) - x(j,1,i))*interforce(x,i,j)*(p(j,1,i) + p(j,1,i+1))/mass(i)
          energy= energy+ p(j,1,i)*siteenergy(p,x,i,j)/mass(i)
       end do
       factors(1)= factors(1)+energy/dble(N)
    end do
    factors(2)=1.0d0

    deallocate(tempx, rp, rk)
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
