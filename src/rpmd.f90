program rpmd
  use potential
  use verletint
  use estimators
  use instantonmod
  use general
  use MKL_VSL_TYPE
  use MKL_VSL

  implicit none
  !general variables
  integer::                        i, j, count,k,jj, nrep
  integer::                        i1,i2,j1,j2,idof1,idof2
  integer::                        idof,ii
  integer::                        time1, time2,irate, imax
  double precision::               answer,sigmaA, weight
  double precision::               tcf0, totalweight,s, ringpot
  double precision, allocatable::  x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  hesstemp(:,:,:,:), v(:,:),p0(:,:,:),q0(:,:,:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:)
  namelist /MCDATA/ n, beta, NMC, noutput,dt, iprint,imin,tau,&
       nrep, use_fft, thermostat, ndim, natom, xunit,gamma, &
       outputtcf, outputfbar

  !-------------------------
  !Set default system parameters then read in namelist
  iprint=.false.
  n= 100
  beta= 100.0d0
  NMC= (5e6)
  Noutput= 1e5
  dt= 1d-3
  nrep=1
  thermostat=1
  ndim=3
  natom=1
  xunit=1
  use_fft=.false.
  imin=0
  tau=1.0d0
  gamma=1.0d0
  nestim=1
  ring=.true.

  read(5, nml=MCDATA)
  betan= beta/dble(n)
  omegan= 1.0/betan
  call system_clock(time1, irate, imax)
  write(*,*)"Running with parameters (in a.u.):"
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "NMC, Noutput, dt=", NMC, Noutput, dt
  if (thermostat .eq. 1) then
     write(*,*) "tau=", tau
     write(*,*) "Running with Andersen thermostat"
  else if (thermostat .eq. 2) then
     write(*,*) "tau=", tau
     write(*,*) "gamma=", gamma
     write(*,*) "Running with Langevin thermostat"
  end if
  call V_init()
  ndof= ndim*natom
  totdof= ndof*n

  allocate(mass(natom))
  open(15,file="mass.dat")
  do i=1,natom
     read(15,*) mass(i)
  end do
  close(15)
  !-------------------------
  !-------------------------
  !Read in the transition state and work out normal
  allocate(transition(ndim, natom), normalvec(ndim,natom))
  open(20, file="transition.dat")
  do i=1, natom
     read(20, *) (transition(j,i), j=1, ndim)
  end do
  close(20)
  if (xunit .eq. 2) transition(:,:)= transition(:,:)/0.529177d0


  allocate(hess(ndof,ndof), hesstemp(ndim, natom, ndim, natom))
  call massweightedhess(transition, hesstemp)
  do i1=1, natom
     do j1=1, ndim
        do i2= 1, natom
           do j2= 1,ndim
              idof1= (j1-1)*ndim + i1
              idof2= (j2-1)*ndim + i2
              hess(idof1,idof2)= hesstemp(j1,i1,j2,i2)
           end do
        end do
     end do
  end do
  lwork=1+ 6*ndof + 2*ndof**2
  liwork= 3 + 5*ndof
  info=0
  allocate(work(lwork), iwork(liwork))
  allocate(transfreqs(ndof))
  call dsyevd('V', 'U', ndof, hess, ndof, transfreqs,work,lwork,iwork,&
       liwork,info)
  deallocate(work, iwork)
  ! write(*,*) hess
  do i=1, natom
     do j=1, ndim
        idof= calcidof(i,j)
        normalvec(j,i)= hess(1,idof)
        if (hess(1,1) .gt. 0.0) then
           normalvec(j,i)= hess(1,idof)
        else
           normalvec(j,i)= -hess(1,idof)
        end if
     end do
  end do
  write(*,*) "normalvec:",normalvec
  !-------------------------
  !-------------------------
  !Start converging TCFs

  call alloc_nm(0)
  call init_nm() !allocates ring polymer stuff
  call init_estimators() !allocates estimator stuff

  allocate(x(n,ndim,natom),p(n,ndim,natom))
  allocate(totaltcf(nestim,ntime), tcf(nestim,ntime))
  if (outputfbar) allocate(p0(n,ndim,natom),q0(n,ndim,natom))
  totaltcf(:,:)=0.0d0
  totalweight=0.0d0

  !--------------------
  !Main loop
  do ii=1, nrep
     tcf(:,:)=0.0d0
     call init_path(x,p, tcf0)
     if (outputfbar) then
        p0(:,:,:)=p(:,:,:)
        q0(:,:,:)=x(:,:,:)
     end if
     call propagator(x,p,tcf)
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf(i,j)*tcf0
        end do
     end do
     !Optional output of f1 and fbar to check stats
     if (outputfbar) then
        allocate(v(ndim,natom))
        do i=1,ndim
           do j=1,natom
              v(i,j)= centroid(x(:,i,j)) - transition(i,j)
           end do
        end do
        s= dot_product(reshape(v,(/ndof/)), reshape(normalvec,(/ndof/)))
        if (s .gt. 0) then
           ringpot= 0.0
           do k=1,n
              ringpot= ringpot+ pot(q0(k,:,:))
           end do
           if (outputfbar) write(200,*) centroid(p0(:,1,1))/mass(1)*exp(-betan*ringpot),&
                p(1,1,1)*exp(-betan*ringpot)/mass(1)

        end if
        deallocate(v)
     end if
  end do

  !------------------------------------
  !Finalize and write out
  totaltcf(:,:)= totaltcf(:,:)/dble(nrep)
  write(*,*) "Reactant partition function:", calcpartition()
  write(*,*) "Final rate:", totaltcf(1,ntime)
  if (outputtcf) then
     open(100, file="timecorrelation.dat")
     do i=1, ntime
        write(100,*) dt*dble(i*noutput), (totaltcf(j,i),j=1,nestim)
     end do
     close(100)
  end if
  call free_nm()
  call finalize_estimators()
  call system_clock(time2, irate, imax)
  write(*,*) "Time taken (s)=",dble(time2-time1)/dble(irate)
  deallocate(mass, hess, hesstemp)
  deallocate(x,p, totaltcf, tcf)
  deallocate(transition, normalvec, transfreqs)

end program rpmd
