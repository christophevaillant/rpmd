module verletint
  use nr_fft
  use potential
  use instantonmod
  use estimators
  use general

  implicit none

  !TODO: sort out public and private

  integer::                         imin, iproc
  integer, allocatable::            ipar(:)
  double precision::                dt, dHdrlimit
  double precision, allocatable::   c1(:,:), c2(:,:)
  double precision::                tau, gamma
  logical::                         iprint
contains

  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory
  subroutine linear_pimd_andersen(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision, intent(in)::  dbdl(:,:),a(:,:),b(:,:)
    double precision, intent(inout)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(out):: dHdr
    double precision::              randno, sigma
    double precision::              totenergy, kin
    double precision, allocatable:: pprop(:), tempp(:), tempv(:)
    integer::              i,j,k,l,count, time1,time2,imax, irate, dofi
    integer (kind=4)::     rkick(1)

    ring=.false.
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
          errcode_poisson = virngpoisson( rmethod_poisson, stream_poisson, 1, rkick(1), dble(Noutput))
          count=0
       end if
       call time_step(xprop, vprop)
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
  !routine taking a step in time using normal mode verlet algorithm
  !from ceriotti et al 2010 paper.
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
  subroutine init_linear(a,b)
    double precision, intent(inout)::    a(:,:),b(:,:)
    integer::             i,j,k,l, dofi

    ring=.false.
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

  end subroutine init_linear
  !-----------------------------------------------------
  !-----------------------------------------------------
  !initialize normal mode routines
  subroutine init_ring()
    integer::             i,j,k,l, dofi

    ring=.true.
    !------------
    !take care of the centroid mode
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
  end subroutine init_ring
  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine linear_pimd_pile(xprop,vprop,a,b,dbdl,dHdr)
    implicit none
    double precision, intent(inout)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(in)::  a(:,:),b(:,:), dbdl(:,:)
    double precision, intent(out):: dHdr
    double precision::              potenergy, springenergy, kinenergy, totenergy
    double precision::              contr,randno, Eold
    integer::              i,j,k,count, time1,time2,imax, irate, skipcount,jj
    integer (kind=4)::     rkick(1)

    ring=.false.
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
             call init_linear_path(a, b, xprop, vprop)
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
  subroutine time_step_ffpile(xprop, vprop)
    implicit none
    double precision, intent(inout)::    xprop(:,:,:), vprop(:,:,:)
    double precision, allocatable::  force(:,:,:), pip(:,:,:)

    allocate(force(n, ndim, natom), pip(n, ndim, natom))
    ! write(*,*) "pot 1:"
    ! write(*,*) "----------"
    call step_v(0.5d0*dt, xprop, vprop, force, .true.)
    ! write(*,*) "thermo:"
    call step_langevin(vprop)
    ! write(*,*) "pot 2:"
    call step_v(0.5d0*dt, xprop, vprop, force,.false.)
    ! write(*,*) "nm:"
    call step_nm(dt,xprop,vprop ,.true.)
    deallocate(force)
    return
  end subroutine time_step_ffpile


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
          ! write(*,*) p(:,i,j)
          do k=1,n
             vprop(k,i,j)= p(k,i,j) !seems silly, but it isn't!
             p(k,i,j)= (c1(j,k)**2)*p(k,i,j) + &
                  sqrt(beadmass(j,k)/betan)*c2(j,k)*sqrt(1.0+c1(j,k)**2)*pprop((dofi-1)*n +k)
             ! write(*,*)k, c1(j,k), c2(j,k)
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
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine ring_rpmd(xprop,vprop,tcf)
    implicit none
    double precision, intent(inout)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(out):: tcf(:,:)
    double precision, allocatable::  x0(:,:,:),q0(:,:,:),q2(:,:,:)
    double precision::              contr,randno, Eold,p0
    integer::              i,j,k,count, time1,time2,imax, irate, jj
    integer (kind=4)::     rkick(1)

    ring=.true.
    allocate(x0(n,ndim,natom))
    x0(:,:,:)=xprop(:,:,:)
    p0= centroid(vprop(:,1,1))
    ! if (use_fft) then
    !    allocate(q0(n,ndim,natom), q2(n,ndim,natom))
    !    call ring_transform_forward_nr(x0(:,1,1),q0(:,1,1))
    !    call ring_transform_forward(xprop(:,1,1),q2(:,1,1))
    !    call ring_transform_backward_nr(q0(:,1,1), x0(:,1,1))
    !    do i=1,n
    !       write(*,*)i,q0(i,1,1),q2(i,1,1)
    !       write(*,*)i,x0(i,1,1), xprop(i,1,1)
    !    end do
    ! stop
    ! end if
    do i=1, NMC, 1
       call time_step(xprop, vprop)
       if (mod(i,Noutput) .eq. 0 .and. i .gt. imin) then !count .ge. Noutput
          call estimator(xprop, vprop, tcf(1,i/Noutput))
          if (iprint) write(*,*) 100*dble(i)/dble(NMC), tcf(1,i/Noutput), centroid(xprop(:,1,1))
       end if
    end do
    deallocate(x0)
    return
  end subroutine ring_rpmd

  !-----------------------------------------------------
  !-----------------------------------------------------
  !main function for propagating a MD trajectory with Langevin thermostat
  subroutine ring_rpmd_pile(xprop,vprop,tcf)
    implicit none
    double precision, intent(inout)::  xprop(:,:,:), vprop(:,:,:)
    double precision, intent(out):: tcf(:,:)
    double precision::              contr,randno, Eold
    integer::              i,j,k,count, time1,time2,imax, irate, jj
    integer (kind=4)::     rkick(1)

    ring=.true.
    allocate(c1(natom,n), c2(natom,n))
    do i=1,n
       do j=1,natom
          ! write(*,*) i, lam(i), mass(j), beadmass(j,i)

          c1(j,i)= exp(-gamma*dt*lam(i)*sqrt(mass(j)/beadmass(j,i)))!
          c2(j,i)= sqrt(1.0d0- c1(j,i)**2)
       end do
    end do
    do i=1, NMC, 1
       ! write(*,*)"---------"
       ! write(*,*) i, xprop
       call time_step_ffpile(xprop, vprop)
       ! write(*,*) xprop
       if (mod(i,Noutput) .eq. 0 .and. i .gt. imin) then !count .ge. Noutput
          call estimator(xprop, vprop, tcf(1,i/Noutput))
          if (iprint) write(*,*) 100*dble(i)/dble(NMC), tcf(1,i/Noutput), centroid(xprop(:,1,1))
       end if
    end do
    deallocate(c1, c2)
    return
  end subroutine ring_rpmd_pile
  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
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
    ! write(*,*) "loc at", imin(1), dists(imin(1))
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
             ! write(*,*) iproc, k,xtilde(1,i,j), xtilde(imin(1), i,j), x(k,i,j)
          end do
       end do
    end do
    deallocate(splinepath, dists)
    ! x(:,:,:)= xtilde(:,:,:)
    ! open(45+iproc)
    ! do i=1,n
    !    write(45+iproc,*) natom
    !    write(45+iproc,*) "Energy of minimum",i
    !    do j=1, natom
    !       write(45+iproc,*)  label(j), (x(i,k,j)*0.529177d0, k=1,ndim)
    !    end do
    ! end do
    ! close(45+iproc)

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
