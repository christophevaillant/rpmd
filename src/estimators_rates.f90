module estimators
  use general
  use instantonmod

  implicit none

  public:: init_estimators, estimator, finalize_estimators, transition, normalvec,&
       transfreqs, hess, init_path

  private

  logical, allocatable:: whichestim(:)
  double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:)
  double precision, allocatable:: lambda(:,:), hess(:,:)

contains
  
  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof

    if (outputfbar) open(200,file="fbar.dat")

    ntime= NMC/Noutput
    allocate(whichestim(nestim))
    ! allocate(lambda(ndof,n))
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
    return
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
    double precision,allocatable:: v(:,:)
    double precision:: s
    integer:: i,j

    allocate(v(ndim,natom))
    do i=1,ndim
       do j=1,natom
          v(i,j)= centroid(x(:,i,j)) - transition(i,j)
       end do
    end do
    s= dot_product(reshape(v,(/ndof/)), reshape(normalvec,(/ndof/)))
    if (s .gt. 0) then
       tcfval=1.0d0
    else
       tcfval=0.0d0
    end if

    deallocate(v)
    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors)
    double precision, intent(out):: x(:,:,:), p(:,:,:),factors
    double precision, allocatable::   vel(:), tempp(:), pos(:), tempx(:)
    double precision, allocatable::  rk(:,:), pk(:), rp(:)
    double precision::                stdev, potvals, ringpot, potdiff, poscent, ringpols
    integer::                         i,j,k, idof, imin(1), inn

    allocate(tempp(n), vel(n),pos(n), rp(n),tempx(n))
    allocate(rk(n,ndof),pk(n))
    potvals=0.0d0
    stdev= 1.0d0/sqrt(betan)
    do idof=1, ndof
          !---------------------------
          !generate random numbers for chain
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,rp,0.0d0,stdev)!ringpolymer
          !---------------------------
          !scale them to have the right standard deviation
          ringpols= 0.0d0
          do k=1, n
             if (k.eq. 1) then
                rp(k)= 0.0d0
             else
                rp(k)= rp(k)/(lam(k))
             end if
          end do
          !---------------------------
          !transform to non-ring polymer normal mode coords
          if (ring) then
             if (use_fft) then
                call ring_transform_backward_nr(rp, rk(:,idof))
             else
                call ring_transform_backward(rp, rk(:,idof))
             end if
          else
             call linear_transform_backward(rp, rk(:,idof), idof)
          end if
          !---------------------------
          !produce random positions in DoFs orthogonal to unstable mode
          if (idof .gt. 1) then
             errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pos,0.0d0,stdev)
             do k=1,n
                rk(k,idof)= rk(k,idof) +pos(k)/sqrt(abs(transfreqs(idof)))
                potvals= potvals+(0.5d0*mass(1)*transfreqs(idof)*pos(k)**2  + V0)/dble(N)
             end do
          end if
    end do
    !---------------------------
    !transform to cartesian coordinates    
    do k=1,n
       x(k,:,:)= reshape(matmul(transpose(hess),rk(k,:)),(/ndim,natom/))
    end do
    do i=1,ndim
       do j=1,natom
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pk(:),0.0d0,stdev)!momenta
          x(:,i,j)= x(:,i,j/sqrt(mass(j)))!mass comes from not having included in ringpolymer and normal mode generation
          p(:,i,j)= p(:,i,j)*sqrt(mass(j))
          !pin centroid to barrier:
          x(:,i,j)= x(:,i,j) - centroid(x(:,i,j)) + transition(i,j)
       end do
    end do
    !calculate real potential for each bead
    ringpot= 0.0d0
    do k=1,n
       ! write(*,*) k, x(k,:,:),pot(x(k,:,:))
       ringpot= ringpot+ pot(x(k,:,:))
    end do
    potdiff= (ringpot- potvals) !potvals is 0 for 1d
    ! write(*,*) ringpot, potvals, potdiff,exp(-betan*potdiff)
    factors=exp(-betan*potdiff)*centroid(p(:,1,1))/mass(1)

    deallocate(vel, tempp,pos, tempx)

    return
  end subroutine init_path

end module estimators
