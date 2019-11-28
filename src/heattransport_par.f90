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
  include 'mpif.h'
  !general variables
  integer::                        i, j, count,k,jj, nrep
  integer::                        i1,i2,j1,j2,idof1,idof2
  integer::                        idof,ii
  integer::                        time1, time2,irate, imax
  double precision::               answer,sigmaA, weight, massin, averagex, testgauss(1)
  double precision::               totalweight,s, ringpot, norm, T, Tleft, Tright
  double precision, allocatable::  x(:,:,:), p(:,:,:), totaltcf(:,:)
  double precision, allocatable::  v(:,:),p0(:,:,:),q0(:,:,:),tcf0(:)
  !lapack variables
  character::                      jobz, range, uplo
  double precision::               vl, vu, abstol, deltaT
  integer::                        nout, ldz, lwork, liwork, info
  integer,allocatable::            isuppz(:), iwork(:)
  double precision, allocatable::  work(:), z(:,:)
  logical::                        latticemass
  !MPI variables
  integer::                        ierr, nproc, ncalcs
  integer (kind=4)::               ierrmkl
  integer, dimension(MPI_STATUS_SIZE) :: rstatus
  integer, allocatable::           calcscounts(:), ptrs(:)
  double precision, allocatable::  alltcfs(:,:,:), allweights(:)

  namelist /MCDATA/ n, beta, NMC, noutput,dt, iprint,imin,tau,&
       nrep, use_fft, thermostat, ndim, natom, xunit,gamma, &
       outputtcf, latticemass, deltaT, convection

  !initialize MPI
  nproc=0
  iproc=0
  ncalcs=0
  ierr=0
  call MPI_Init(ierr)
  !get the processor ID
  call MPI_Comm_rank ( MPI_COMM_WORLD, iproc, ierr )
  !get the total number of processors
  call MPI_Comm_size ( MPI_COMM_WORLD, nproc, ierr )

  !-------------------------
  !Set default system parameters then read in namelist
  iprint=.false.
  n= 100
  beta= 100.0d0
  NMC= (5e6)
  Noutput= 1e5
  dt= 1d-3
  nrep=1
  thermostat=0
  ndim=3
  natom=1
  xunit=1
  use_fft=.false.
  imin=0
  tau=1.0d0
  gamma=1.0d0
  nestim=2
  ring=.true.
  fixedends=.false.
  latticemass=.false.
  deltaT=0.0d0
  outputtcf=.true.
  convection=.false.

  !-----------------------------------------
  !Read in parameters for processor 0 and distribute to all processors
  if (iproc .eq. 0) then
     read(5, nml=MCDATA)
     betan= beta/dble(n)
     omegan= 1.0d0/betan
     call system_clock(time1, irate, imax)
     write(*,*)"Running with parameters (in a.u.):"
     write(*,*) "beta, betan, n=", beta, betan, n
     write(*,*) "ndim, natom=", ndim, natom
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
  end if

  call MPI_Bcast(nrep, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(N, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(NMC, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(thermostat, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(imin, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(ndim, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(natom, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(convection, 1, MPI_LOGICAL, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(use_fft, 1, MPI_LOGICAL, 0,MPI_COMM_WORLD, ierr)

  call V_init(iproc)

  ndof= ndim*natom
  totdof= ndof*n
  Noutput=1 !force noutput to enforce baths at every step
  allocate(mass(natom))
  allocate(transition(ndim, natom), normalvec(ndim,natom))
  allocate(hess(ndof,ndof), transfreqs(ndof))

  if (iproc .eq.0) then
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
     !Read in the transition state and work out normal
     open(20, file="transition.dat")
     do i=1, natom
        read(20, *) (transition(j,i), j=1, ndim)
     end do
     close(20)
     if (xunit .eq. 2) transition(:,:)= transition(:,:)/0.529177d0

     call findhess(transition, hess)
     call findnormal(hess, transfreqs,normalvec)
     write(*,*) "Harmonic frequencies:",transfreqs
  end if

  !-------------------------------------
  !broadcast to all processors
  call MPI_Bcast(beta, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mass, natom, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(spacing, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(harm, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(interharm, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(gammaheat, 1, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(lattice, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(transition, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(transfreqs, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(normalvec, ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  call MPI_Bcast(hess, ndof*ndof, MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD, ierr)
  ! call MPI_Bcast(label, natom, MPI_CHARACTER, 0,MPI_COMM_WORLD, ierr) !might need this later

  !-----------------------------
  !collection of variables to be set after broadcast
  betan= beta/dble(n)
  omegan= 1.0d0/betan
  V0= pot(transition)

  T= 1.0d0/(3.167d-6*beta)
  Tleft= T- 0.5d0*deltaT
  Tright= T+ 0.5d0*deltaT

  betaleft=1.0d0/(3.167d-6*Tleft)
  betaright=1.0d0/(3.167d-6*Tright)

  ncalcs= int(nrep/nproc)
  if (iproc .eq. 0) ncalcs= ncalcs + mod(nrep,nproc)
  write(*,*) iproc, "Energy zero=", V0, "Calcs:", ncalcs
  !-------------------------
  !-------------------------
  !Start converging TCFs

  call alloc_nm(iproc)
  call init_nm() !allocates ring polymer stuff
  call init_estimators() !allocates estimator stuff

  allocate(x(n,ndim,natom),p(n,ndim,natom))
  allocate(totaltcf(nestim,ntime), tcf(nestim,ntime), tcf0(nestim))
  totaltcf(:,:)=0.0d0
  totalweight=0.0d0
  !--------------------
  !Main loop

  do ii=1, ncalcs
     tcf(:,:)=0.0d0
     weight=0.0d0
     x(:,:,:)=0.0d0
     p(:,:,:)=0.0d0
     call init_path(x,p, tcf0, weight)
     totalweight=totalweight + weight
     if (iprint .and. iproc .eq. 0) write(*,*) ii,tcf0(1), weight, totalweight/dble(ii)
     call propagator(x,p,tcf)
     do i=1, nestim
        do j=1,ntime
           totaltcf(i,j)= totaltcf(i,j) + tcf(i,j)*tcf0(i)*weight
        end do
     end do
  end do
  totaltcf(:,:)= totaltcf(:,:)/totalweight
  totalweight=totalweight/dble(ncalcs)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  !-------------------------------------------------
  !Collect the tcfs and weights from all processors
  allocate(alltcfs(nestim, ntime,nproc), allweights(nproc), calcscounts(nproc), ptrs(nproc))
  alltcfs(:,:,:)=0.0d0
  allweights=0.0d0
  calcscounts(:)= nestim*ntime
  do i=1, nproc
     ptrs(i)= (i-1)*nestim*ntime
  end do
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Gatherv(totaltcf, nestim*ntime, MPI_DOUBLE_PRECISION, alltcfs, &
       calcscounts,ptrs, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Gather(totalweight, 1, MPI_DOUBLE_PRECISION, allweights, 1, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !------------------------------------
  !Finalize and write out
  if (iproc .eq. 0) then
     write(*,*) ntime, NMC, noutput
     totaltcf(:,:)=0.0d0
     totalweight=0.0d0
     do ii=1, nproc
        do i=1, nestim
           do j=1, ntime
              totaltcf(i,j)= totaltcf(i,j) + alltcfs(i,j,ii)/dble(nproc)
           end do
        end do
        totalweight= totalweight+ allweights(ii)/dble(nproc)
     end do
     write(*,*) "totalweight:", totalweight
     write(*,*) allweights(:)
     open(100, file="timecorrelation.dat")
     do i=1, ntime
        write(100,*) dt*dble(i*noutput), (totaltcf(j,i),j=1,nestim)
     end do
     close(100)
     call system_clock(time2, irate, imax)
     write(*,*) "Time taken (s)=",dble(time2-time1)/dble(irate)
  end if
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call free_nm()
  call finalize_estimators()
  deallocate(mass)
  deallocate(x,p, totaltcf, tcf)
  deallocate(hess, normalvec)

end program rpmd
