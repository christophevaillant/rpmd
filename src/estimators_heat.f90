module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,transfreqs, init_path, whichestim,&
       betanright, betanleft, width

  private

  logical, allocatable:: whichestim(:)
  double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:)
  double precision, allocatable:: lambda(:,:)
  double precision:: betanleft, betanright, width

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
    double precision:: s, energy, prob, stdev
    integer:: i,j,k

    allocate(newp(1))
    tcfval(1)=0.0d0
    allocate(grad(N,1,natom))
    call UMprime(x, grad)
    do i=1, natom
       energy=0.0d0
       do j=1, N

       end do
       tcfval(1)= tcfval(1) + energy/dble(N)
    end do
    deallocate(newp,grad)

    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors)
    double precision, intent(out):: x(:,:,:), p(:,:,:),factors
    double precision, allocatable:: tempx(:),rk(:,:),rp(:), grad(:,:,:)
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
    do i=1,ndim
       do j=1,natom
          idof= calcidof(i,j)
          stdevharm= 1.0d0/sqrt(harm*mass(j))
          stdev= 1.0d0/sqrt(betan)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,1,tempx,lattice(i,j),stdevharm)!sites
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,p(:,i,j),0.0d0,stdev)!momenta
          x(:,i,j)= rk(:,idof) + tempx(1)
          p(:,i,j)= p(:,i,j)*sqrt(mass(j))
       end do
    end do
    potvals=0.0d0
    do i=1,ndim
       do j=1,natom-1
          do k=1,n
             potvals=potvals + 0.5d0*mass(j)*(interharm*(x(k,i,j)-x(k,i,j+1)))**2
          end do
       end do
    end do
    ringpot= 0.0d0
    allocate(grad(N,1,natom))
    do k=1,n
       ringpot= ringpot+ pot(x(k,:,:))
    end do
    potdiff= (ringpot- potvals)/dble(natom)
    ! write(*,*) ringpot, potvals, potdiff
    factors=0.0d0
    call UMprime(x, grad)
    do i=1,natom
       energy=0.0d0
       do j=1,n

       end do
       factors=factors+ energy/dble(n)
    end do
    factors=factors*exp(-betan*potdiff)
    deallocate(tempx, rp, rk,grad)
    return
  end subroutine init_path

end module estimators
