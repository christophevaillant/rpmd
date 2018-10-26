program rpmd
  use potential
  use verletint
  use estimators
  use instantonmod
  use MKL_VSL_TYPE
  use MKL_VSL

  implicit none
  !general variables
  integer::                        i, j, count,k,jj, nrep
  integer::                        i1,i2,j1,j2,idof1,idof2
  integer::                        thermostat,idof,ii
  integer::                        time1, time2,irate, imax
  double precision::               answer,sigmaA, weight, totalweight
  double precision, allocatable::  x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  hesstemp(:,:,:,:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:)
  namelist /MCDATA/ n, beta, NMC, noutput,dt, iprint,imin,tau,&
       nrep, use_mkl, thermostat, ndim, natom, xunit,gamma
       
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
  use_mkl=.false.
  imin=0
  tau=1.0d0
  gamma=1.0d0
  nestim=1

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
  
  !-------------------------
  !-------------------------
  !Read in the transition state and work out normal
  allocate(transition(ndim, natom), normalvec(ndim,natom))
  open(20, file="transition.dat")
  do i=1, natom
     read(20, *) (transition(j,i), j=1, ndim)
  end do
  if (xunit .eq. 2) transition(:,:)= transition(:,:)/0.529177d0
     
  allocate(hess(ndof,ndof), hesstemp(ndim, natom, ndim, natom))
  call Vdoubleprime(transition, hesstemp)
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
  do i=1, natom
     do j=1, ndim
        idof= (j1-1)*ndim + i1
        normalvec(j,i)= hess(1,idof)
     end do
  end do
              
  !-------------------------
  !-------------------------
  !Start converging TCFs
  
  call alloc_nm(0)
  call init_ring() !allocates ring polymer stuff
  call init_estimators() !allocates estimator stuff

  allocate(x(n,ndim,natom),p(n,ndim,natom))
  allocate(totaltcf(nestim,ntime), tcf(nestim,ntime))
  totaltcf(:,:)=0.0d0
  totalweight=0.0d0
  do ii=1, nrep
     call init_path(x,p, weight)
     totalweight= totalweight+weight
     call ring_rpmd_pile(x,p,tcf)
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf(i,j)*weight
        end do
     end do
  end do
  
  totaltcf(:,:)= totaltcf(:,:)/totalweight
  open(100, file="timecorrelation.dat")
  do i=1, ntime
     write(100,*) dt*dble(i*noutput)/dble(nmc), (tcf(j,i),j=1,nestim)
  end do
  deallocate(x,p, totaltcf)
  deallocate(transition, normalvec)
end program rpmd
