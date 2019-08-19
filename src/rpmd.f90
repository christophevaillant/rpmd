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
  double precision::               answer,sigmaA, weight, averagex, qr
  double precision::               tcf0, totalweight,s, ringpot
  double precision, allocatable::  x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  v(:,:),p0(:,:,:),q0(:,:,:)
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


  allocate(hess(ndof,ndof), transfreqs(ndof))
  call findhess(transition, hess)
  call findnormal(transition, hess, transfreqs,normalvec)
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
  averagex=0.0d0
  !--------------------
  !Main loop
  qr=calcpartition()
  do ii=1, nrep
     tcf(:,:)=0.0d0
     call init_path(x,p, tcf0, weight)
     if(centroid(p(:,1,1)).gt.0)  averagex=averagex+tcf0*weight
     if (outputfbar) then
        p0(:,:,:)=p(:,:,:)
        q0(:,:,:)=x(:,:,:)
     end if
     totalweight=totalweight+weight
     call propagator(x,p,tcf)
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf(i,j)*tcf0*weight
        end do
     end do
     if (iprint) write(*,*) ii,centroid(x(:,1,1)),centroid(p(:,1,1)), tcf0, weight, totalweight/dble(ii), totaltcf(1,ntime)/dble(ii)
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
  write(*,*) "average x test:", averagex/totalweight
  !------------------------------------
  !Finalize and write out
  totaltcf(:,:)= totaltcf(:,:)/dble(nrep) !partition function cancels in 1D
  write(*,*) "Beta, betan=", beta, betan
  write(*,*) "Reactant partition function:", qr
  write(*,*) "Final rate:", totaltcf(1,ntime)
  write(*,*) "Final rate (ms^-1):", 2.187730d6*totaltcf(1,ntime)
  write(*,*) "Total weight:", totalweight
  write(*,*) "QTST estimate:", totaltcf(1,1)
  if (outputtcf) then
     open(100, file="timecorrelation.dat")
     do i=1, ntime
        write(100,*) dt*dble(i*noutput), (totaltcf(j,i),j=1,nestim)
     end do
     close(100)
  end if
  if (ndof.eq.1) then
     write(*,*) "classical TST=", exp(-beta*V0)/sqrt(2.0d0*pi*beta*mass(1))
     write(*,*) "classical correction factor=",totaltcf(1,ntime)*sqrt(2.0d0*pi*beta*mass(1))*exp(beta*V0)
  end if
  call free_nm()
  call finalize_estimators()
  call system_clock(time2, irate, imax)
  write(*,*) "Time taken (s)=",dble(time2-time1)/dble(irate)
  deallocate(mass, hess)
  deallocate(x,p, totaltcf, tcf)
  deallocate(transition, normalvec, transfreqs)

end program rpmd
