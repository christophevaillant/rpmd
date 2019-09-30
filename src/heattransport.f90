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
  double precision::               answer,sigmaA, weight, massin
  double precision::               tcf0, totalweight,s, ringpot
  double precision, allocatable::  x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  v(:,:),p0(:,:,:),q0(:,:,:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol, deltabeta
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:)
  logical::                        latticemass
  namelist /MCDATA/ n, beta, NMC, noutput,dt, iprint,imin,tau,&
       nrep, use_fft, thermostat, ndim, natom, xunit,gamma, &
       outputtcf, latticemass, deltabeta, width

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
  fixedends=.false.
  latticemass=.false.
  deltabeta=1.0d0
  outputtcf=.true.
  width=1.0d-1

  read(5, nml=MCDATA)
  betan= beta/dble(n)
  omegan= 1.0/betan
  call system_clock(time1, irate, imax)
  write(*,*)"Running with parameters (in a.u.):"
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "NMC, Noutput, dt=", NMC, Noutput, dt
  write(*,*) "Nrep=", nrep
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

  betanleft=(beta-deltabeta)/dble(n)
  betanright=(beta+deltabeta)/dble(n)

  Noutput=1 !force noutput to enforce baths at every step

  allocate(mass(natom))
  open(15,file="mass.dat")
  if (latticemass) then
     read(15,*) massin
     do i=1,natom
        mass(i)=massin
     end do
  else
     do i=1,natom
        read(15,*) mass(i)
     end do
  end if
  close(15)

  !-------------------------
  !-------------------------
  !Start converging TCFs

  call alloc_nm(0)
  call init_nm() !allocates ring polymer stuff
  call init_estimators() !allocates estimator stuff

  allocate(x(n,ndim,natom),p(n,ndim,natom))
  allocate(totaltcf(nestim,ntime), tcf(nestim,ntime))
  totaltcf(:,:)=0.0d0
  totalweight=0.0d0

  !--------------------
  !Main loop
  do ii=1, nrep
     tcf(:,:)=0.0d0
     call init_path(x,p, tcf0)
     if (iprint) write(*,*) ii, x(1,1,1), p(1,1,1),tcf0
     call propagator(x,p,tcf)
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf(i,j)*tcf0
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
  deallocate(mass)
  deallocate(x,p, totaltcf, tcf)

end program rpmd