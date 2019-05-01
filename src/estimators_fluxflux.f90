module estimators
  use general
  use instantonmod

  implicit none

  public:: init_estimators, estimator, finalize_estimators, surface1, surface2, &
       freqs1, freqs2, hess1,hess2, init_path, normalvec1, normalvec2, width

  private

  logical, allocatable::          whichestim(:)
  double precision, allocatable:: surface1(:,:), normalvec1(:,:), freqs1(:)
  double precision, allocatable:: surface2(:,:), normalvec2(:,:), freqs2(:)
  double precision, allocatable:: lambda(:,:), hess1(:,:), hess2(:,:)
  double precision::              width

contains

  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof
    !TODO: estimatormod- implement list of estimators

    if (outputfbar) open(200,file="fbar.dat")

    ntime= NMC/Noutput
    allocate(whichestim(nestim))
    allocate(lambda(ndof,n))
    ! do idof=1, ndof
    !    do k=1, n,2
    !       kk=(k-1)/2
    !       if (idof .eq. 1) then !skip unstable mode
    !          if (kk.eq.0) then
    !             lambda(idof, k)= (2.0*omegan*sin(pi*kk/N))**2
    !          else
    !             lambda(idof, k-1)= (2.0*omegan*sin(pi*kk/N))**2
    !             lambda(idof, k)= (2.0*omegan*sin(pi*kk/N))**2
    !          end if
    !       else
    !          lambda(idof, k)= ((4.0*omegan**2 + 2.0*transfreqs(idof)) &
    !               - 4.0*omegan**2*cos(2.0*pi*kk/n))
    !          if (k.gt.1) &
    !               lambda(idof, k-1)= ((4.0*omegan**2 + 2.0*transfreqs(idof)) &
    !               - 4.0*omegan**2*cos(2.0*pi*kk/n))
    !       end if
    !       if (mod(n,2).eq.0) then
    !          ! kk=n/2
    !          ! if (idof .eq. 1 .and. ndof .eq. 1) then
    !          ! lambda(idof, n)= (2.0*omegan*sin(pi*kk/N))**2
    !          ! else
    !          lambda(idof, n)= ((4.0*omegan**2 + 2.0*transfreqs(idof)/mass(1)) &
    !               - 4.0*omegan**2*cos(2.0*pi*kk/n))
    !       end if
    !    end do
    ! end do
    ! write(*,*) lambda
    ! write(*,*) transfreqs
  end subroutine init_estimators

  subroutine finalize_estimators()
    implicit none
    if (outputfbar) close(200)
    deallocate(whichestim)
    deallocate(lambda)
  end subroutine finalize_estimators

  subroutine estimator(x,p,tcfval)
    implicit none
    double precision, intent(in):: x(:,:,:), p(:,:,:)
    double precision, intent(out):: tcfval(:)
    double precision,allocatable:: v(:,:,:)
    double precision:: s, prob
    integer:: i,j

    tcfval(:)=0.0d0
    do i=1, n
       prob= exp(-0.5d0*(x(i,1,1)-surface2(1,1))**2/width**2)/sqrt(2.0d0*pi*width**2)
       tcfval(1)= tcfval(1)+ p(i,1,1)*prob/mass(1)
       tcfval(2)= tcfval(2)+ prob
    end do
    tcfval(1)= tcfval(1)/dble(n)
    tcfval(2)= tcfval(2)/dble(n)
    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors)
    double precision, intent(out):: x(:,:,:), p(:,:,:),factors(:)
    double precision, allocatable::   vel(:), tempp(:), pos(:), tempx(:)
    double precision, allocatable::  rk(:,:), pk(:,:), rp(:)
    double precision::                stdev, potvals, ringpot, potdiff, poscent
    integer::                         i,j,k, idof, imin(1), inn

    allocate(tempp(n), vel(n),pos(n), rp(n),tempx(n))
    allocate(rk(n,ndof),pk(n,ndof))
    potvals=0.0d0
    stdev= 1.0d0/sqrt(betan)
    do idof=1, ndof
          !---------------------------
          !generate random numbers for chain
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,rp,0.0d0,stdev)!ringpolymer
          ! errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,vel,0.0d0,stdev)!momenta
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
          !calculate contribution from the harmonic potential (ignore
          ! the contribution from the bead springs)
          ! if (idof .gt. 1) then
          !    do k=1,n
          !       potvals= potvals+0.5d0*transfreqs(idof)*rp(k)**2
          !    end do
          ! end if
          !---------------------------
          !transform to non-ring polymer normal mode coords
          if (ring) then
             if (use_fft) then
                call ring_transform_backward_nr(rp, rk(:,idof))
             else
                call ring_transform_backward(rp, rk(:,idof))
             end if
             ! call ring_transform_backward(vel, pk(:,idof))
          else
             call linear_transform_backward(rp, rk(:,idof), idof)
             ! call linear_transform_backward(vel, pk(:,idof), 0)
          end if
          if (idof .gt. 1) then
             errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pos,0.0d0,stdev)
             do k=1,n
                rk(k,idof)= rk(k,idof) +pos(k)/sqrt(freqs1(idof))
             end do
          end if
          rk(:,idof)= rk(:,idof) - rk(1,idof) + surface1(1,1)
    end do
    !---------------------------
    !transform to cartesian coordinates    
    ringpot= 0.0
    do k=1,n
       x(k,:,:)= reshape(matmul(hess1,rk(k,:)),(/ndim,natom/))
       ! p(k,:,:)= reshape(matmul(transpose(hess),pk(k,:)),(/ndim,natom/))
    end do
    do i=1,ndim
       do j=1,natom
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,p(:,i,j),0.0d0,stdev)!momenta
          x(:,i,j)= (x(:,i,j)/sqrt(mass(j))) + surface1(i,j)
          p(:,i,j)= p(:,i,j)*sqrt(mass(j))
       end do
    end do
    do k=1,n
       ringpot= ringpot+ pot(x(k,:,:))
    end do
    potdiff= (ringpot- potvals)
    factors(1)=0.5*exp(-betan*potdiff)*p(1,1,1)/mass(1) !half because half the chain has to be between surface 1 and surface2, so not all beads are equivalent
    factors(2)=0.0d0
    do i=1,ndim
       do j=1,natom
          factors(2)= factors(2) + (normalvec1(i,j)*normalvec2(i,j))**2
       end do
    end do
    factors(2)= exp(-betan*potdiff)*sqrt(factors(2))
    deallocate(vel, tempp,pos, tempx)
    return
  end subroutine init_path

end module estimators
