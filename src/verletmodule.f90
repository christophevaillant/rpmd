module verletint
  use potential
  use instantonmod
  use estimators
  use general

  implicit none

  !TODO: sort out public and private

  public::   propagator, thermostat, dt, imin, iproc, tau, gamma, iprint,&
       init_nm

  private

  integer::                         imin, iproc, thermostat
  integer, allocatable::            ipar(:)
  double precision::                dt, dHdrlimit
  double precision, allocatable::   c1(:,:), c2(:,:)
  double precision::                tau, gamma
  logical::                         iprint
contains

  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Andersen thermostat
  subroutine propagator(xprop,vprop,tcf,a,b,dbdl)
    implicit none
    double precision, intent(inout)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(out):: tcf(:,:)
    double precision::              sigma
    double precision, allocatable:: pprop(:), tempp(:), tempv(:)
    integer::              dofi, i,j,k,l,count
    integer (kind=4)::     rkick(1)
    double precision, intent(in), optional::  dbdl(:,:),a(:,:),b(:,:)

    !----------------------------------------
    !Initialize  and allocate
    allocate(pprop(n), tempp(n),tempv(n))

    count=0
    tcf(:,:)=0.0d0
    if (thermostat .eq. 1) then
       errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
    else if (thermostat .eq. 2) then
       allocate(c1(natom,n), c2(natom,n))
       do i=1,n
          do j=1,natom
             c1(j,i)= exp(-gamma*dt*lam(i)*sqrt(mass(j)/beadmass(j,i)))!
             c2(j,i)= sqrt(1.0d0- c1(j,i)**2)
          end do
       end do
    end if

    do i=1, NMC, 1
       !----------------------------------------
       !Deal with thermostats
       if (thermostat.eq.1) then !andersen thermostat
          count=count+1
          if (count .ge. rkick(1)) then
             count=0
             ! if (iprint) write(*,*) 100*dble(i)/dble(NMC),tcf(1,i/Noutput)
             do j=1,ndim
                do k=1,natom
                   dofi=calcidof(j,k)
                   sigma=sqrt(1.0d0/betan)
                   errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pprop,0.0d0,sigma)
                   do l=1,n
                      tempp(l)= pprop(l)*sqrt(beadmass(k,l))
                   end do
                   if (.not. use_fft) then
                      if (ring) then
                         call ring_transform_backward(tempp, tempv)
                      else
                         call linear_transform_backward(tempp, tempv, 0)
                      end if
                   else
                      if (ring) then
                         call ring_transform_backward_nr(tempp, tempv)
                      else
                         call linear_transform_backward_nr(tempp, tempv,0)
                      end if
                   end if
                   do l=1,n
                      vprop(l, j, k) = tempv(l)
                   end do
                end do
             end do
             errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, &
                  rkick(1), dble(Noutput))
             count=0
          end if
       end if
       !----------------------------------------
       !Take a step
       if (thermostat .eq. 2) then !langevin thermostat
          call time_step_ffpile(xprop, vprop)
       else
          call time_step(xprop, vprop) !no thermostat
       end if
       !----------------------------------------
       !Calculate the estimator for this step
       if (i .gt. imin) then
          call estimator(xprop, vprop, tcf(:,i/Noutput))
       end if
    end do
    deallocate(pprop, tempv, tempp)
    if (thermostat .eq. 2) deallocate(c1,c2)
    !TODO: add estimator_finalize
    return
  end subroutine propagator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step(xprop, vprop)
    implicit none
    double precision, intent(inout)::    xprop(:,:,:), vprop(:,:,:)
    double precision, allocatable::  force(:,:,:), pip(:,:,:)

    allocate(force(n, ndim, natom), pip(n, ndim, natom))
    call step_v(dt, xprop, vprop, force, .true.)
    call step_nm(dt,xprop,vprop ,.true.)
    deallocate(force)
    return
  end subroutine time_step
  !-----------------------------------------------------
  !-----------------------------------------------------
  !initialize normal mode routines
  subroutine init_nm(a,b)
    double precision, intent(inout), optional::    a(:,:),b(:,:)
    integer::             i,j,k,l, dofi

    if (.not. ring) then
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
    else
       k=0
       lam(1)= 0.0d0
       do j=1, natom
          beadmass(j,1)=mass(j)
       end do
       do l=1,n
          transmatrix(1,l)= 1.0d0/sqrt(dble(n)) !k=0
       end do
       !do all the other modes in steps of two
       do i=2, n-2,2
          k= (i/2)
          lam(i)= 2.0d0*sin(dble(k)*pi/dble(n))/betan
          lam(i+1)=2.0d0*sin(dble(k)*pi/dble(n))/betan
          do j=1, natom
             !TODO: implement possible mass scalings
             beadmass(j,i)=mass(j)!*(lam(i)*tau)**2
             beadmass(j,i+1)=mass(j)!*(lam(i)*tau)**2
          end do
          do l=1,n
             transmatrix(i,l)= sqrt(2.0d0/dble(n))*&
                  cos(2.0d0*dble((l-1)*k)*pi/dble(n))
             transmatrix(i+1,l)= sqrt(2.0d0/dble(n))*&
                  sin(2.0d0*dble((l-1)*k)*pi/dble(n))
             if (transmatrix(i,l) .ne. transmatrix(i,l)) then
                write(*,*)"Nan!"
                stop
             end if
          end do
       end do
       if (mod(n,2).eq.0) then
          i=n
          k=(i/2)
          lam(i)=2.0d0/betan
          do j=1,natom
             beadmass(j,i)=mass(j)!*(lam(i)*tau)**2
          end do
          do l=0,n-1
             transmatrix(i,l+1)= ((-1.0d0)**l)/sqrt(dble(n))
          end do
       end if

    end if
  end subroutine init_nm

  !-----------------------------------------------------
  !-----------------------------------------------------
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
  subroutine time_step_ffpile(xprop, vprop)
    implicit none
    double precision, intent(inout)::    xprop(:,:,:), vprop(:,:,:)
    double precision, allocatable::  force(:,:,:), pip(:,:,:)

    allocate(force(n, ndim, natom), pip(n, ndim, natom))
    call step_v(0.5d0*dt, xprop, vprop, force, .true.)
    call step_langevin(vprop)
    call step_v(0.5d0*dt, xprop, vprop, force,.false.)
    call step_nm(dt,xprop,vprop ,.true.)
    deallocate(force)
    return
  end subroutine time_step_ffpile

  !-----------------------------------------------------
  !-----------------------------------------------------
  !Normal mode steps
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
          if (.not. use_fft) then
             if (ring) then
                call ring_transform_forward(p(:,i,j), pip(:,i,j))
                call ring_transform_forward(x(:,i,j), q(:,i,j))
             else
                call linear_transform_forward(p(:,i,j), pip(:,i,j), 0)
                call linear_transform_forward(x(:,i,j), q(:,i,j), dofi)
             end if
          else
             if (ring) then
                call ring_transform_forward_nr(p(:,i,j), pip(:,i,j))
                call ring_transform_forward_nr(x(:,i,j), q(:,i,j))
             else
                call linear_transform_forward_nr(p(:,i,j), pip(:,i,j), 0)
                call linear_transform_forward_nr(x(:,i,j), q(:,i,j), dofi)
             end if
          end if
       end do
    end do
    do i=1, n
       do j=1,ndim
          do k=1,natom
             if (i.eq.1 .and. ring) then
                newpi(i,j,k)=pip(i,j,k)
                q(i,j,k)= q(i,j,k) + &
                     pip(i,j,k)*time/beadmass(k,i)
             else
                omegak= sqrt(mass(k)/beadmass(k,i))*lam(i)
                newpi(i,j,k)= pip(i,j,k)*cos(time*omegak) - &
                     q(i,j,k)*omegak*beadmass(k,i)*sin(omegak*time)
                q(i,j,k)= q(i,j,k)*cos(time*omegak) + &
                     pip(i,j,k)*sin(omegak*time)/(omegak*beadmass(k,i))
             end if
             if (newpi(i,j,k) .ne. newpi(i,j,k)) then
                write(*,*) "NaN in 1st NM propagation"
                stop
             end if
          end do
       end do
    end do

    if (transform) then
       do i=1,ndim
          do j=1,natom
             dofi=calcidof(i,j)
             if (.not. use_fft) then
                if (ring) then
                   call ring_transform_backward(newpi(:,i,j), p(:,i,j))
                   call ring_transform_backward(q(:,i,j), x(:,i,j))
                else
                   call linear_transform_backward(newpi(:,i,j), p(:,i,j), 0)
                   call linear_transform_backward(q(:,i,j), x(:,i,j),dofi)
                end if
             else
                if (ring) then
                   call ring_transform_backward_nr(newpi(:,i,j), p(:,i,j))
                   call ring_transform_backward_nr(q(:,i,j), x(:,i,j))
                else
                   call linear_transform_backward_nr(newpi(:,i,j), p(:,i,j),0)
                   call linear_transform_backward_nr(q(:,i,j), x(:,i,j),dofi)
                end if
             end if
          end do
       end do
    else
       x(:,:,:)= q(:,:,:)
       p(:,:,:)= newpi(:,:,:)
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
             p(i,j,k)= p(i,j,k) - force(i,j,k)*time
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

  subroutine step_langevin(vprop)
    implicit none
    double precision, intent(inout):: vprop(:,:,:)
    double precision, allocatable:: p(:,:,:), pprop(:)
    integer:: i,j,k, dofi

    allocate(p(n,ndim,natom), pprop(n))
    errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n*ndof,pprop,0.0d0,1.0d0)
    do i=1,ndim
       do j=1,natom
          dofi= calcidof(i,j)
          if (.not. use_fft) then
             if (ring) then
                call ring_transform_forward(vprop(:,i,j), p(:,i,j))
             else
                call linear_transform_forward(vprop(:,i,j), p(:,i,j), 0)
             end if
          else
             if (ring) then
                call ring_transform_forward_nr(vprop(:,i,j), p(:,i,j))
             else
                call linear_transform_forward_nr(vprop(:,i,j), p(:,i,j), 0)
             end if
          end if
          do k=1,n
             vprop(k,i,j)= p(k,i,j) !seems silly, but it isn't!
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pprop((dofi-1)*n +k)
          end do
       end do
    end do
    do k=1,n
       do j=1,natom
          p(k,:,j)= norm2(p(k,:,j))*vprop(k,:,j)/norm2(vprop(k,:,j)) !see, not so silly!
       end do
    end do
    do i=1,ndim
       do j=1,natom
          dofi= (j-1)*ndim +i
          if (.not. use_fft) then
             if (ring) then
                call ring_transform_backward(p(:,i,j),vprop(:,i,j))
             else
                call linear_transform_backward(p(:,i,j),vprop(:,i,j), 0)
             end if
          else
             if (ring) then
                call ring_transform_backward_nr(p(:,i,j),vprop(:,i,j))
             else
                call linear_transform_backward_nr(p(:,i,j),vprop(:,i,j), 0)
             end if
          end if
       end do
    end do

  end subroutine step_langevin

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  !TODO: move init_linear_path to a new estimators module, where relevant
  subroutine init_linear_path(start, end, x, p)
    double precision::  start(:,:), end(:,:), x(:,:,:), p(:,:,:)
    double precision, allocatable:: vel(:), tempp(:), dists(:)
    double precision, allocatable:: splinepath(:)
    double precision::  stdev
    integer::           i,j,k, dofi, imin(1)

    allocate(vel(n),tempp(n), dists(n))

    do i=1,n
       dists(i)= eucliddist(x(i,:,:), end(:,:))
    end do
    imin= minloc(dists)
    do i=1, imin(1)
       dists(i)= dble(i-1)/dble(imin(1)-1)
    end do
    allocate(splinepath(imin(1)))
    do i=1,ndim
       do j=1,natom
          splinepath(:)=0.0d0
          call spline(dists(1:imin(1)), x(1:imin(1),i,j), 1.0d31, 1.0d31, splinepath(:))
          do k=1,n
             x(k,i,j)= splint(dists(1:imin(1)), x(1:imin(1),i,j), splinepath(:), dble(k-1)/dble(n-1))
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
  end subroutine init_linear_path

end module verletint
