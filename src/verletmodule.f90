include 'mkl_vsl.fi'
module verletint
  use potential
  use instantonmod
  use MKL_VSL_TYPE
  use MKL_VSL
  implicit none

  integer::                         NMC, Noutput, imin, iproc
  integer, allocatable::            ipar(:)
  double precision::                dt, dHdrlimit
  double precision,allocatable::    transmatrix(:,:),beadvec(:,:)
  double precision,allocatable::    xtilde(:,:,:),beadmass(:,:), lam(:)
  double precision, allocatable::   c1(:,:), c2(:,:)
  double precision::                tau, gamma
  logical::                         iprint, use_mkl
  integer::                         errcode_poisson, rmethod_poisson
  integer::                         brng_poisson, seed_poisson, stat
  type (vsl_stream_state)::         stream_poisson,stream_normal
  integer::                         errcode_normal, rmethod_normal, brng_normal
  integer::                         seed_normal, rk

contains
  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(start, end, x, p)
    double precision, intent(in)::    start(:,:), end(:,:)
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision, allocatable::   vel(:), tempp(:), dists(:)
    double precision, allocatable::   splinepath(:)
    double precision::                stdev
    integer::                         i,j,k, dofi, imin(1)

    allocate(vel(n),tempp(n), dists(n))

    do i=1,n
       dists(i)= eucliddist(xtilde(i,:,:), end(:,:))
    end do
    imin= minloc(dists)
    do i=1, imin(1)
       dists(i)= dble(i-1)/dble(imin(1)-1)
    end do
    allocate(splinepath(imin(1)))
    do i=1,ndim
       do j=1,natom
          splinepath(:)=0.0d0
          call spline(dists(1:imin(1)), xtilde(1:imin(1),i,j), 1.0d31, 1.0d31, splinepath(:))
          do k=1,n
             x(k,i,j)= splint(dists(1:imin(1)), xtilde(1:imin(1),i,j), splinepath(:), dble(k-1)/dble(n-1))
          end do
       end do
    end do
    deallocate(splinepath, dists)

    do i=1,ndim
       do k=1,natom
          dofi=(k-1)*ndim +i
          stdev=sqrt(1.0d0/betan)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,vel,0.0d0,stdev)
          do j=1,n
             vel(j)= vel(j)*sqrt(beadmass(k,j))
          end do
          call linear_transform_backward(vel, tempp, 0)
          do j=1,n
             p(j,i,k)= tempp(j)
          end do
       end do
    end do

    deallocate(vel, tempp)
    return
  end subroutine init_path

  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory
  subroutine linear_pimd_andersen(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision, intent(in)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(in)::  dbdl,a(:,:),b(:,:)
    double precision, intent(out):: dHdr
    double precision::              randno, sigma
    double precision::              totenergy, kin,dbdl(:,:)
    double precision, allocatable:: pprop(:), tempp(:), tempv(:)
    integer::              i,j,k,l,count, time1,time2,imax, irate, dofi
    integer (kind=4)::     rkick(1)

    allocate(pprop(n), tempp(n),tempv(n))
    count=0
    dHdr=0.0d0
    kin=0.0d0
    errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
    do i=1, NMC, 1
       count=count+1
       if (count .ge. rkick(1)) then
          count=0
          if (iprint) write(*,*) 100*dble(i)/dble(NMC),dHdr/(betan**2*dble(i-imin))
          do j=1,ndim
             do k=1,natom
                dofi=(k-1)*ndim +j
                sigma=sqrt(1.0d0/betan)
                errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pprop,0.0d0,sigma)
                do l=1,n
                   tempp(l)= pprop(l)*sqrt(beadmass(k,l))
                end do
                if (.not. use_mkl) then
                   call linear_transform_backward(tempp, tempv, 0)
                else
                   call linear_transform_backward_nr(tempp, tempv,0)
                end if
                do l=1,n
                   vprop(l, j, k) = tempv(l)
                end do
             end do
          end do
          errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
          count=0
       end if
       call time_step_test(xprop, vprop)
       if (i .gt. imin) then
       do j=1,ndim
          do k=1,natom
             dHdr= dHdr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
          end do
       end do
    end if
    end do
    dHdr=dHdr/dble(NMC-imin)
    deallocate(pprop, tempv, tempp)
    return
  end subroutine linear_pimd_andersen
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode forward transformation for linear polymer
  subroutine linear_transform_forward(xprop, qprop, bead)
    implicit none
    double precision, intent(in)::  xprop(:)
    double precision, intent(out):: qprop(:)
    integer, intent(in)::           bead
    integer::                 i,j
    qprop(:)=xprop(:)
    call dsymv('U', n, 1.0d0,transmatrix, n, xprop,1,0.0d0, qprop,1)
    if (bead .gt. 0) then
       do i=1,n
          qprop(i)= qprop(i) - beadvec(i,bead)
       end do
    end if
    return
  end subroutine linear_transform_forward
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode backward transformation for linear polymer
  subroutine linear_transform_backward(qprop, xprop,bead)
    implicit none
    double precision, intent(out):: xprop(:)
    double precision, intent(in)::  qprop(:)
    integer, intent(in)::           bead
    integer::                 i,j
    if (bead .gt. 0) then
       do i=1,n
          qprop(i)= qprop(i) + beadvec(i,bead)
       end do
    end if
    call dsymv('U', n, 1.0d0, transmatrix, n, qprop,1,0.0d0, xprop,1)
    if (bead .gt.0) then
       do i=1,n
          qprop(i)= qprop(i) - beadvec(i,bead)
       end do
    end if
    return
  end subroutine linear_transform_backward

  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step(x, v)
    implicit none
    double precision, intent(inout):: x(:,:,:), v(:,:,:)
    double precision, allocatable::   newv(:,:,:), newx(:,:,:),force(:,:)
    double precision, allocatable::   q(:,:,:), p(:,:,:)
    integer::             i,j,k,dofi
    double precision::    omegak

    allocate(newv(n,ndim,natom), newx(n,ndim,natom),force(ndim,natom))
    allocate(q(n,ndim,natom),p(n,ndim,natom))
    !-------------
    !step 1
    do i=1,n
       call Vprime(x(i,:,:),force)
       do j=1,ndim
          do k=1,natom
             newv(i,j,k)= v(i,j,k) - force(j,k)*dt
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 2
    do i=1,ndim
       do j=1,natom
          if (.not. use_mkl) then
             dofi= (j-1)*ndim+i
             call linear_transform_forward(newv(:,i,j), p(:,i,j), 0)
             call linear_transform_forward(x(:,i,j), q(:,i,j), dofi)
          else
             call linear_transform_forward_nr(newv(:,i,j), p(:,i,j), 0)
             call linear_transform_forward_nr(x(:,i,j), q(:,i,j), dofi)
          end if
       end do
    end do
    !-------------
    !step 3
    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newv(i,j,k)= p(i,j,k)*cos(dt*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*dt)
             newx(i,j,k)= q(i,j,k)*cos(dt*omegak) + &
                  p(i,j,k)*sin(omegak*dt)/(omegak*beadmass(k,i))
             if (newv(i,j,k) .ne. newv(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do
    !-------------
    !step 4
    do i=1,ndim
       do j=1,natom
          dofi=(j-1)*ndim + i
          if (.not. use_mkl) then
             call linear_transform_backward(newv(:,i,j), v(:,i,j), 0)
             call linear_transform_backward(newx(:,i,j), x(:,i,j),dofi)
          else
             call linear_transform_backward_nr(newv(:,i,j), v(:,i,j),0)
             call linear_transform_backward_nr(newx(:,i,j), x(:,i,j),dofi)
          end if
       end do
    end do

    deallocate(newv, newx,force)
    deallocate(q,p)
    return
  end subroutine time_step_test
  !-----------------------------------------------------
  !-----------------------------------------------------
  !initialize normal mode routines
  subroutine init_nm(a,b)
    double precision, intent(inout)::    a(:,:),b(:,:)
    integer::             i,j,k,l, dofi

    do i=1, n
       lam(i)= 2.0d0*sin(dble(i)*pi/dble(2*n +2))/betan
       do j=1, natom
          beadmass(j,i)=mass(j)*(lam(i)*tau)**2
       end do
       do l=i,n
          transmatrix(i,l)= sin(dble(i*l)*pi/dble(n+1))*sqrt(2.0d0/dble(n+1))
          transmatrix(l,i)= transmatrix(i,l)
          if (transmatrix(i,l) .ne. transmatrix(i,l)) then
             write(*,*)"Nan!"
             stop
          end if
       end do
       do j=1,ndim
          do k=1,natom
             dofi= (k-1)*ndim +j
             beadvec(i, dofi)= a(j,k)*sin(dble(i)*pi/dble(n+1)) &
                  + b(j,k)*sin(dble(n*i)*pi/dble(n+1))
             beadvec(i,dofi)= beadvec(i,dofi)*sqrt(2.0d0/dble(n+1))
             beadvec(i,dofi)= beadvec(i,dofi)/(lam(i)*betan)**2
          end do
       end do
    end do

  end subroutine init_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !free normal mode arrays

  subroutine free_nm()
    deallocate(transmatrix,beadmass,beadvec,lam)
  end subroutine free_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !allocate normal mode arrays

  subroutine alloc_nm(iproc)
    implicit none
    integer::    itime, irate, imax
    integer, intent(in):: iproc

    call system_clock(itime,irate,imax)
    seed_normal= mod(itime+5*iproc,1000)
    call system_clock(itime,irate,imax)
    seed_poisson= mod(itime+5*iproc,1000)
    brng_normal = VSL_BRNG_MT19937
    brng_poisson = VSL_BRNG_MT19937
    rmethod_normal = VSL_RNG_METHOD_GAUSSIAN_ICDF
    rmethod_poisson = VSL_RNG_METHOD_POISSON_POISNORM
    errcode_normal = vslnewstream( stream_normal, brng_normal, seed_normal )
    errcode_poisson = vslnewstream( stream_poisson, brng_poisson, seed_poisson )
    
    write(*,*)"proc", iproc, "running with seeds:", seed_normal, seed_poisson

  end subroutine alloc_nm
  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine linear_pimd_pile(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision, intent(in)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(in)::  a(:,:),b(:,:), dbdl(:,:)
    double precision, intent(out):: dHdr
    double precision::              potenergy, springenergy, kinenergy, totenergy
    double precision::              contr,randno, Eold
    integer::              i,j,k,count, time1,time2,imax, irate, skipcount,jj
    integer (kind=4)::     rkick(1)

    allocate(c1(natom,n), c2(natom,n))
    do i=1,n
       do j=1,natom
          c1(j,i)= exp(-gamma*dt*lam(i)*sqrt(mass(j)/beadmass(j,i)))!
          c2(j,i)= sqrt(1.0d0- c1(j,i)**2)
       end do
    end do
    count=0
    dHdr=0.0d0
    skipcount=0
    do i=1, NMC, 1
       count=count+1
       if (count .ge. Noutput .and. i .gt. imin) then
          count=0
          if (iprint) write(*,*) 100*dble(i)/dble(NMC-imin), dHdr/(betan**2*dble(i-imin))
       end if
       call time_step_ffpile(xprop, vprop)
       if (i.gt.imin) then
          contr=0.0d0
          do j=1,ndim
             do k=1,natom
                contr=contr+mass(k)*(-xprop(n,j,k))*dbdl(j,k)
             end do
          end do
          if (abs(contr) .lt. dHdrlimit .or. dHdrlimit .lt. 0.0) then
             dHdr= dHdr+contr
          else
             write(*,*) "Over limit", contr, ", reinitialize path"
             call init_path(a, b, xprop, vprop)
          end if
       end if
    end do
    dHdr= dHdr/dble(NMC-imin-skipcount)
    ! if (skipcount .gt. 0) write(*,*) "Skipped", skipcount, "points"
    deallocate(c1, c2)
    return
  end subroutine linear_pimd_pile

  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_ffpile(xprop, vprop, force)
    implicit none
    double precision::    xprop(:,:,:), vprop(:,:,:), omegak
    double precision, allocatable::  newv(:,:,:), newx(:,:,:)
    double precision, allocatable::  q(:,:,:), p(:,:,:)
    double precision, allocatable::  pprop(:),force(:,:,:)
    integer::             i,j,k,dofi

    allocate(force(n, ndim, natom))

    call step_v(0.5d0*dt, xprop, vprop, force, .true.)
    call step_langevin(vprop)
    call step_v(0.5d0*dt, xprop, vprop, force,.false.)
    call step_nm(dt,xprop,vprop ,.true.)

    deallocate(force)
    return
  end subroutine time_step_ffpile

  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode forward transformation for linear polymer
  subroutine linear_transform_forward_nr(xprop, qprop, bead)
    implicit none
    double precision::        xprop(:), qprop(:)
    double precision, allocatable:: qprop_nr(:)
    integer::                 i,j,bead
    allocate(qprop_nr(1:n+1))
    qprop_nr(1)=0.0d0
    qprop_nr(2:n+1)=xprop(1:n)
    call sinft(qprop_nr)
    qprop(1:n)= qprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
    if (bead.gt.0) qprop(:)= qprop(:) - beadvec(:,bead)
    deallocate(qprop_nr)
    return
  end subroutine linear_transform_forward_nr
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode backward transformation for linear polymer
  subroutine linear_transform_backward_nr(qprop, xprop,bead)
    implicit none
    double precision::        xprop(:), qprop(:)
    double precision, allocatable:: xprop_nr(:)
    integer::                 i,j,bead
    allocate(xprop_nr(1:n+1))
    if (bead .gt. 0) then
       xprop_nr(2:n+1)=qprop(1:n) + beadvec(1:n, bead)
    else 
       xprop_nr(2:n+1)=qprop(1:n)
    end if
    call sinft(xprop_nr)
    xprop(1:n)= xprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
    deallocate(xprop_nr)
    return
  end subroutine linear_transform_backward_nr


  !-----------------------------------------------------
  !-----------------------------------------------------
  subroutine step_nm(time, x, p, transform)
    implicit none
    integer::          i, j, k, dofi
    double precision:: omegak, newpi(n, ndim, natom)
    double precision:: q(n, ndim, natom), pip(n, ndim, natom)
    double precision:: x(:,:,:), p(:,:,:) , time
    logical::          transform

       do i=1,ndim
          do j=1,natom
             dofi= (j-1)*ndim+i
             if (.not. use_mkl) then
                call linear_transform_forward(p(:,i,j), pip(:,i,j), 0)
                call linear_transform_forward(x(:,i,j), q(:,i,j), dofi)
             else
                call linear_transform_forward_nr(p(:,i,j), pip(:,i,j), 0)
                call linear_transform_forward_nr(x(:,i,j), q(:,i,j), dofi)
             end if
          end do
       end do

    do i=1, n
       do j=1,ndim
          do k=1,natom
             omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
             newpi(i,j,k)= pip(i,j,k)*cos(time*omegak) - &
                  q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*time)
             q(i,j,k)= q(i,j,k)*cos(time*omegak) + &
                  pip(i,j,k)*sin(omegak*time)/(omegak*beadmass(k,i))
             if (newpi(i,j,k) .ne. newpi(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do
    pip(:,:,:)=newpi(:,:,:)
    
    if (transform) then
       do i=1,ndim
          do j=1,natom
             dofi=(j-1)*ndim + i
             if (.not. use_mkl) then
                call linear_transform_backward(pip(:,i,j), p(:,i,j), 0)
                call linear_transform_backward(q(:,i,j), x(:,i,j),dofi)
             else
                call linear_transform_backward_nr(pip(:,i,j), p(:,i,j),0)
                call linear_transform_backward_nr(q(:,i,j), x(:,i,j),dofi)
             end if
          end do
       end do
    else
       x(:,:,:)= q(:,:,:)
       p(:,:,:)= pip(:,:,:)
    end if

       
  end subroutine step_nm
  !-----------------------------------------------------
  !-----------------------------------------------------

  subroutine step_v(time,x,p, force, recalculate)
    implicit none
    double precision, intent(in):: force(:, :,:), time
    double precision, intent(inout):: p(:,:,:), x(:,:,:)
    integer::         i,j,k
    logical, intent(in):: recalculate

    do i=1,n
       if (recalculate) call Vprime(x(i,:,:),force(i,:,:))
       do j=1,ndim
          do k=1,natom
             p(i,j,k)= p(i,j,k) - force(j,k)*time
             if (p(i,j,k) .ne. p(i,j,k)) then
                write(*,*) "NaN in pot propagation"
                stop
             end if
          end do
       end do
    end do

  end subroutine step_v
  !-----------------------------------------------------
  !-----------------------------------------------------

  subroutine step_langevin(p)
    implicit none
    double precision, intent(inout):: v(:,:,:)
    double precision, allocatable:: p(:,:,:)
    integer:: i,j,k, dofi

    !-------------
    !Thermostat step
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pprop,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim+i
          if (.not. use_mkl) then
             call linear_transform_forward(v(:,i,j), p(:,i,j), 0)
          else
             call linear_transform_forward_nr(v(:,i,j), p(:,i,j), 0)
          end if
          do k=1,n
             v(k,i,j)= p(k,i,j) !seems silly, but it isn't!
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pprop((dofi-1)*n +k)
          end do
       end do
    end do
    do k=1,n
       do j=1,natom
          p(k,:,j)= norm2(p(k,:,j))*v(k,:,j)/norm2(v(k,:,j))
       end do
    end do
    
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim +i
          if (.not. use_mkl) then
             call linear_transform_backward(p(:,i,j),vprop(:,i,j), 0)
          else
             call linear_transform_backward_nr(p(:,i,j),vprop(:,i,j), 0)
          end if
       end do
    end do

  end subroutine step_langevin
end module verletint
