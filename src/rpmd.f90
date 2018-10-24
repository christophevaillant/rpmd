program rpmd
  use potential
  use verletint
  use estimators
  use MKL_VSL_TYPE
  use MKL_VSL

  implicit none
  !general variables
  integer::                        i, j, count,k,jj, nrep
  integer::                        thermostat,dofi
  double precision::               dHdr, answer,sigmaA
  double precision, allocatable::  origin(:), wellinit(:,:), vgrad(:,:)
  double precision, allocatable::  x(:,:,:), vel(:), tempv(:), tempp(:), initpath(:,:)
  double precision, allocatable::  p(:,:,:), xi(:), dbdxi(:,:,:), Vpath(:)
  double precision, allocatable::  path(:,:,:), lampath(:), splinepath(:), y(:,:,:)
  integer::                        ndofrb, dummy
  character::                      dummylabel, dummystr(28)
  !gauss-legendre variables
  integer::                        nintegral,ii, time1, time2, imax, irate
  double precision, allocatable::  weights(:), xint(:,:,:), integrand(:), sigma(:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  diag(:), work(:), z(:,:)
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
  nintegral= 5
  nrep=1
  thermostat=1
  ndim=3
  natom=1
  xunit=1
  use_mkl=.false.
  imin=0
  tau=1.0d0
  gamma=1.0d0

  read(5, nml=MCDATA)
  betan= beta/dble(n+1)
  call system_clock(time1, irate, imax)
  write(*,*)"Running with parameters (in a.u.):"
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "NMC, Noutput, dt=", NMC, Noutput, dt
  write(*,*) "with integration points", nintegral
  write(*,*) "tau=", tau
  write(*,*) "gamma=", gamma
  if (thermostat .eq. 1) then
     write(*,*) "Running with Andersen thermostat"
  else if (thermostat .eq. 2) then
     write(*,*) "Running with Langevin thermostat"
  end if
  call V_init()

  !-------------------------
  !-------------------------
  !Start loop over integration points
  allocate(integrand(nintegral), sigma(nintegral))
  integrand(:)=0.0d0
  call random_seed()
  allocate(x(n,ndim,natom), vel(n),p(n,ndim,natom))
  allocate(tempp(n), tempv(n))
  allocate(transmatrix(n,n),beadmass(natom,n),beadvec(n,ndof), lam(n))
  call alloc_nm(1)
  sigma(ii)=0.0d0
  integrand(ii)=0.0d0
  do jj=1, nrep
     call init_nm(path(1,:,:),xint(ii,:,:))
     call init_path(path(1,:,:), xint(ii,:,:), x, p)
     if (thermostat .eq. 1) then
        call propagate_pimd_nm(x,p, path(1,:,:),xint(ii,:,:),dbdxi(ii,:,:),dHdr)
     else if (thermostat .eq. 2) then
        call propagate_pimd_pile(x,p, path(1,:,:),xint(ii,:,:), dbdxi(ii,:,:),dHdr)
     else if (thermostat .eq. 3) then
        call propagate_pimd_higher(x,p, path(1,:,:),xint(ii,:,:),dbdxi(ii,:,:),dHdr)
     else
        write(*,*) "Incorrect thermostat option."
        stop
     end if
     integrand(ii)=integrand(ii) + dHdr/(dble(nrep)*betan**2)
     sigma(ii)= sigma(ii)+ (dHdr/betan**2)**2
  end do
  sigma(ii)= sigma(ii)/dble(nrep)
  sigma(ii)= sigma(ii) - integrand(ii)**2
  ! write(*,*)"-----------------"
  write(*,*) ii, xi(ii), integrand(ii), sqrt(sigma(ii))
  ! write(*,*)"-----------------"
  call free_nm()
  deallocate(x, vel,p, tempv, tempp)
  answer= 0.0d0
  sigmaA=0.0d0

  do i=1, nintegral
     if (integrand(i) .eq. integrand(i)) &
          answer= answer+ weights(i)*integrand(i)
     sigmaA= sigmaA + sigma(i)*weights(i)**2
  end do

  finalI= exp(-answer*betan)
  write(*,*) "beta, betan, n=", beta, betan, n
  write(*,*) "final Delta A=", answer , "+/-", sqrt(sigmaA)/sqrt(dble(nrep))
  write(*,*) "q/q0=", finalI, betan*sqrt(sigmaA)*finalI/sqrt(dble(nrep))
  write(*,*) "q0/q=", 1.0d0/finalI, betan*sqrt(sigmaA)/(finalI*sqrt(dble(nrep)))
  call system_clock(time2, irate, imax)
  write(*,*) "Ran in", dble(time2-time1)/dble(irate), "s"
  deallocate(path, xtilde)
end program rpmd
