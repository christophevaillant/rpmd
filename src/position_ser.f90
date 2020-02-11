program rpmd
  use potential
  use verletint
  use estimators
  use instantonmod
  use general
  use estimators
  use MKL_VSL_TYPE
  use MKL_VSL

  implicit none
  !general variables
  integer::                        i, j, count,k,jj, nrep
  integer::                        i1,i2,j1,j2,idof1,idof2
  integer::                        idof,ii
  integer::                        time1, time2,irate, imax
  double precision::               answer,sigmaA, weight, massin,Z,imweight
  double precision::               totalweight,s, ringpot, norm, T, Tleft, Tright
  double precision, allocatable::  x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  v(:,:),p0(:,:,:),q0(:,:,:),tcf0(:),averagex(:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol
  integer::                        nout, ldz, lwork, liwork, info
  namelist /MCDATA/ n, beta, NMC, noutput,dt, iprint,imin,tau,&
       nrep, use_fft, thermostat, ndim, natom, xunit,gamma, &
       outputtcf, nonlinear

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
  ring=.true.
  fixedends=.false.
  outputtcf=.true.
  nonlinear=0

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
  else if (thermostat .eq. 3) then
     write(*,*) "Running with Matsubara propagator"
  end if

  if (nonlinear .eq. 0) then
     write(*,*) "Running without nonlinear sampling."
  else if (nonlinear .eq. 1) then
     write(*,*) "Running with regular nonlinear sampling."
  else if (nonlinear .eq. 2) then
     write(*,*) "Running with averaged nonlinear sampling."
  end if
  ! if (mod(splitbead,n).ne.0) then
  !    write(*,*) "Wrong number of beads for splitting scheme specified:",splitbead, mod(splitbead,n)
  !    stop
  ! end if

  call V_init(0)
  ndof= ndim*natom
  totdof= ndof*n


  allocate(mass(1))
  open(15,file="mass.dat")
  read(15,*) mass(1)
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
  call findnormal(hess, transfreqs,normalvec)
  write(*,*) "Freqs:", transfreqs
  write(*,*) "Normals:", normalvec
  V0= pot(transition)
  write(*,*) "Energy zero=", V0
  !-------------------------
  !-------------------------
  !Start converging TCFs

  call alloc_nm(0)
  call init_nm() !allocates ring polymer stuff
  call init_estimators() !allocates estimator stuff

  allocate(x(n,ndim,natom),p(n,ndim,natom),tcf0(nestim))
  allocate(totaltcf(nestim,ntime), tcf(nestim,ntime), averagex(nestim))
  totaltcf(:,:)=0.0d0
  totalweight=0.0d0
  averagex(:)=0.0d0
  tcf0(:)=0.0d0

  Z=0.0d0
  do i=1,ndof
     Z=Z+0.5d0/sinh(0.5d0*sqrt(transfreqs(i))*beta)
  end do

  !--------------------
  !Main loop

  do ii=1, nrep
     tcf(:,:)=0.0d0
     call init_path(x,p, tcf0, weight)
     do i=1,nestim
        averagex(i)= averagex(i)+ weight*tcf0(i)**2
        tcf(i,1)= tcf0(i)
     end do
     if (iprint .and. mod(ii,Noutput) .eq. 0 ) write(*,*) ii,centroid(p(:,1,1)), weight, totalweight/dble(ii), averagex(:)/totalweight
     call propagator(x,p,tcf)
     totalweight=totalweight + weight
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf0(i)*tcf(i,j)*weight
        end do
     end do
  end do

  !------------------------------------
  !Finalize and write out
  write(*,*) "Zero time::", averagex(:)/totalweight
  write(*,*) "Harmonic partition:", Z
  write(*,*) "Totalweight:", totalweight/dble(nrep)

  if (nonlinear.eq.0) then
     totaltcf(:,:)= totaltcf(:,:)/totalweight
  else
     do i=1, nestim/2
        totaltcf(i,:)= totaltcf(i,:)*averagex(i+3)/(totalweight*totaltcf(i,1))
     end do
  end if

  if (outputtcf) then
     open(100, file="timecorrelation.dat")
     do i=1, ntime
        write(100,*) dt*dble(i), (totaltcf(j,i),j=1,nestim)
     end do
     close(100)
  end if
  call free_nm()
  call finalize_estimators()
  call system_clock(time2, irate, imax)
  write(*,*) "Time taken (s)=",dble(time2-time1)/dble(irate)

  deallocate(mass)
  deallocate(x,p, totaltcf, tcf)
  deallocate(hess, normalvec)

end program rpmd
