program fluxflux
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
  double precision::               totalweight,s, ringpot
  double precision, allocatable::  tcf0(:),x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  hesstemp(:,:,:,:), v(:,:),p0(:,:,:),q0(:,:,:)
  !lapack variables
  namelist /MCDATA/ n, beta, NMC, noutput,dt, iprint,imin,tau,&
       nrep, use_fft, thermostat, ndim, natom, xunit,gamma, &
       outputtcf, outputfbar, width

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
  nestim=2
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
  allocate(surface1(ndim, natom), normalvec1(ndim,natom))
  allocate(surface2(ndim, natom), normalvec2(ndim,natom))
  open(20, file="surface1.dat")
  do i=1, natom
     read(20, *) (surface1(j,i), j=1, ndim)
  end do
  close(20)
  open(21, file="surface2.dat")
  do i=1, natom
     read(21, *) (surface2(j,i), j=1, ndim)
  end do
  close(21)
  if (xunit .eq. 2) then
     surface1(:,:)= surface1(:,:)/0.529177d0
     surface2(:,:)= surface2(:,:)/0.529177d0
  end if

  allocate(hess1(ndof,ndof),hess2(ndof,ndof))
  allocate(freqs1(ndof),freqs2(ndof))
  call findhess(surface1,hess1)
  call findhess(surface2,hess2)
  call findnormal(surface1, hess1, freqs1, normalvec1)
  call findnormal(surface2, hess2, freqs2, normalvec2)

  write(*,*) "normalvec1:",normalvec1
  write(*,*) "normalvec2:",normalvec2
  !-------------------------
  !-------------------------
  !Start converging TCFs

  call alloc_nm(0)
  call init_nm() !allocates ring polymer stuff
  call init_estimators() !allocates estimator stuff

  allocate(x(n,ndim,natom),p(n,ndim,natom))
  allocate(totaltcf(nestim,ntime), tcf(nestim,ntime), tcf0(nestim))
  totaltcf(:,:)=0.0d0
  totalweight=0.0d0

  !--------------------
  !Main loop
  do ii=1, nrep
     tcf(:,:)=0.0d0
     call init_path(x,p, tcf0)
     if (iprint) write(*,*) ii, tcf0(1), tcf0(2)
     call propagator(x,p,tcf)
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf(i,j)*tcf0(i)
        end do
     end do
  end do

  !------------------------------------
  !Finalize and write out
  totaltcf(:,:)= totaltcf(:,:)/dble(nrep)
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
  deallocate(mass, hess1, hess2)
  deallocate(x,p, totaltcf, tcf)
  deallocate(surface1, surface2)
  deallocate(normalvec1, normalvec2, freqs1, freqs2)

end program fluxflux
