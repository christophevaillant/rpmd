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

    !----------------------------
    !TODO: Not sure if we'll need this for the heat transport. If not, delete.
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
    ! write(*,*) lambda
    ! write(*,*) transfreqs
  end subroutine init_estimators

  subroutine finalize_estimators()
    implicit none
    deallocate(whichestim)
  end subroutine finalize_estimators

  subroutine estimator(x,p,tcfval)
    implicit none
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision, intent(out):: tcfval(:)
    double precision,allocatable:: v(:,:), newp(:)
    double precision:: s, energy, prob, stdev
    integer:: i,j,k

    !TODO: multidimensional will need to work out which side of bath each particle is
    !Keeping this bit of code for future reference, could use to work out box
    ! allocate(v(ndim,natom))
    ! do i=1,ndim
    !    do j=1,natom
    !       v(i,j)= centroid(x(:,i,j)) - transition(i,j)
    !    end do
    ! end do
    ! s= dot_product(reshape(v,(/ndof/)), reshape(normalvec,(/ndof/)))
    ! deallocate(v)

    allocate(newp(1)) !TODO: multidimensional generalization to ndim
    tcfval(1)=0.0d0
    do i=1, natom
       energy=0.0d0
       do j=1, N
          prob= exp(-0.5d0*(x(j,1,i)-lattice(1,i+1))**2/width**2)/sqrt(2.0d0*pi*width**2)
          energy= energy + prob*0.5d0*p(j,1,i)**3/mass(i)**2
          !check the atom is still in the box!
          if (x(j,1,i) .lt. lattice(1,1)) then
             ! write(*,*) "left bounce", x(j,1,i), lattice(1,1)
             x(j,1,i)= lattice(1,1)
             stdev= 1.0d0/sqrt(betanleft) !TODO: check that this shouldn't just be beta
             errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,1,newp,0.0d0,stdev)
             p(j,1,i)= abs(newp(1))*sqrt(mass(i))
          else if (x(j,1,i) .gt. lattice(1,natom+2)) then
             ! write(*,*) "right bounce", x(j,1,i), lattice(1,natom+2)
             x(j,1,i)= lattice(1,Natom+2)
             stdev= 1.0d0/sqrt(betanright) !TODO: check that this shouldn't just be beta
             errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,1,newp,0.0d0,stdev)
             !TODO: multidimensional generalization like for the directional langevin thermostat
             p(j,1,i)= -abs(newp(1))*sqrt(mass(i))
          end if
       end do
       tcfval(1)= tcfval(1) + energy/dble(N)
    end do
    

    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors)
    double precision, intent(out):: x(:,:,:), p(:,:,:),factors
    double precision, allocatable:: tempx(:),rk(:,:),rp(:)
    double precision::              stdev, potvals, ringpot, potdiff
    double precision::              prob, stdevharm, energy
    integer::                       i,j,k, idof, imin(1), inn

    allocate(rp(n),tempx(n))
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
          stdevharm= 1.0d0/sqrt(harm*mass(j))
          stdev= 1.0d0/sqrt(betan)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,tempx,lattice(i,j+1),stdevharm)!sites
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,p(:,i,j),0.0d0,stdev)!momenta
          x(:,i,j)= x(:,i,j) + tempx(:)
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
    do k=1,n
       ringpot= ringpot+ pot(x(k,:,:))
    end do
    potdiff= (ringpot- potvals)
    write(*,*)potdiff, ringpot, potvals
    factors=0.0d0
    do i=1,natom
       energy=0.0d0
       do j=1,n
          prob= exp(-0.5d0*(x(j,1,i)-lattice(1,i+1))**2/width**2)/sqrt(2.0d0*pi*width**2)
          energy= energy + prob*0.5d0*p(j,1,i)**3/mass(i)**2
       end do
       factors=factors+ energy/dble(n)
    end do
    factors=factors*min(exp(-betan*potdiff),1.0d0)
    deallocate(tempx, rp, rk)
    return
  end subroutine init_path

end module estimators
