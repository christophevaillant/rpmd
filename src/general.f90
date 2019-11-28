  include 'mkl_vsl.fi'
  module general
    use MKL_VSL_TYPE
    use MKL_VSL
    use MKL_DFTI
    ! use nr_fft

    implicit none

    integer::                         NMC, Noutput,nestim, ntime
    integer::                         errcode_poisson, rmethod_poisson
    integer::                         brng_poisson, seed_poisson
    type (vsl_stream_state)::         stream_poisson,stream_normal
    type(dfti_descriptor), pointer :: fft_handle
    integer::                         errcode_normal, rmethod_normal, brng_normal
    integer::                         seed_normal, lensav, lenwork
    integer::                         n, ndim, ndof, natom, xunit, totdof
    double precision, allocatable::   tcf(:,:), work(:)
    double precision,allocatable::    transmatrix(:,:),beadvec(:,:)
    double precision,allocatable::    beadmass(:,:), lam(:), wsave(:)
    logical::                         ring, use_fft
    double precision, parameter::    pi=3.14159265358979d0
    double precision::               beta, betan, UMtilde, omegan, dt
    double precision, allocatable::  well1(:,:), well2(:,:), mass(:)
    double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:), hess(:,:)

    character, allocatable::         label(:)
    logical::                        fixedends, outputtcf, outputfbar

    private:: wsave, lensav, work
  contains

    !-----------------------------------------------------
    !-----------------------------------------------------
    !allocate normal mode arrays

    subroutine alloc_nm(iproc)
      implicit none
      integer::    itime, irate, imax, ier
      integer, intent(in):: iproc

      allocate(transmatrix(n,n), beadmass(natom,n), lam(n))
      if (use_fft) then
         ier= DftiCreateDescriptor(fft_handle, dfti_double,&
              dfti_real, 1, n)
         ier= dftisetvalue(fft_handle, dfti_packed_format, dfti_pack_format)
         ier= dftisetvalue(fft_handle, dfti_backward_scale, 1.00d0/sqrt(dble(n)))
         ier= dftisetvalue(fft_handle, dfti_forward_scale, 1.00d0/sqrt(dble(n)))
         ier= DftiCommitDescriptor(fft_handle)
      end if
      if (.not. ring) allocate(beadvec(n,ndof))
      call system_clock(itime,irate,imax)
      seed_normal= mod(itime+5*iproc,1000)
      call system_clock(itime,irate,imax)
      seed_poisson= mod(itime+5*iproc,1000)
      brng_normal = VSL_BRNG_MT19937
      brng_poisson = VSL_BRNG_MT19937
      rmethod_normal = VSL_RNG_METHOD_GAUSSIAN_ICDF
      rmethod_poisson = VSL_RNG_METHOD_POISSON_POISNORM
      errcode_normal = vslnewstream( stream_normal, brng_normal, seed_normal )
      errcode_poisson = vslnewstream( stream_poisson, brng_poisson, seed_poisson )

      write(*,*)"proc", iproc, "running with seeds:", seed_normal, seed_poisson

    end subroutine alloc_nm

    !-----------------------------------------------------
    !-----------------------------------------------------
    !free normal mode arrays

    subroutine free_nm()
      implicit none
      integer:: ier
      deallocate(transmatrix,beadmass,lam)
      if (.not. ring) deallocate(beadvec)
      if (use_fft) ier = DftiFreeDescriptor(fft_handle)
    end subroutine free_nm

    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode forward transformation for linear polymer
    subroutine linear_transform_forward(xprop, qprop, bead)
      implicit none
      double precision, intent(in)::  xprop(:)
      double precision, intent(out):: qprop(:)
      integer, intent(in)::           bead
      integer::                 i,j
      qprop(:)=xprop(:)
      call dsymv('U', n, 1.0d0,transmatrix, n, xprop,1,0.0d0, qprop,1)
      if (bead .gt. 0) then
         do i=1,n
            qprop(i)= qprop(i) - beadvec(i,bead)
         end do
      end if
      return
    end subroutine linear_transform_forward
    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode backward transformation for linear polymer
    subroutine linear_transform_backward(qprop, xprop,bead)
      implicit none
      double precision, intent(out):: xprop(:)
      double precision, intent(inout)::  qprop(:)
      integer, intent(in)::           bead
      integer::                 i,j
      if (bead .gt. 0) then
         do i=1,n
            qprop(i)= qprop(i) + beadvec(i,bead)
         end do
      end if
      call dsymv('U', n, 1.0d0, transmatrix, n, qprop,1,0.0d0, xprop,1)
      if (bead .gt.0) then
         do i=1,n
            qprop(i)= qprop(i) - beadvec(i,bead)
         end do
      end if
      return
    end subroutine linear_transform_backward
    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode forward transformation for ring polymer
    subroutine ring_transform_forward(xprop, qprop)
      implicit none
      double precision, intent(in)::  xprop(:)
      double precision, intent(out):: qprop(:)
      integer::                 i,j,k

      call dgemv('N', n,n, 1.0d0,transmatrix, n, xprop,1,0.0d0, qprop,1)

      return
    end subroutine ring_transform_forward
    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode backward transformation for linear polymer
    subroutine ring_transform_backward(qprop, xprop)
      implicit none
      double precision, intent(out):: xprop(:)
      double precision, intent(in)::  qprop(:)
      integer::                 i,j

      call dgemv('T', n, n,1.0d0, transmatrix, n, qprop,1,0.0d0, xprop,1)
      return
    end subroutine ring_transform_backward

    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode forward transformation for linear polymer
    subroutine linear_transform_forward_nr(xprop, qprop, bead)
      implicit none
      double precision::        xprop(:), qprop(:)
      double precision, allocatable:: qprop_nr(:)
      integer::                 i,j,bead
      allocate(qprop_nr(1:n+1))
      qprop_nr(1)=0.0d0
      qprop_nr(2:n+1)=xprop(1:n)
      ! call sinft(qprop_nr)
      qprop(1:n)= qprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
      if (bead.gt.0) qprop(:)= qprop(:) - beadvec(:,bead)
      deallocate(qprop_nr)
      return
    end subroutine linear_transform_forward_nr
    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode backward transformation for linear polymer
    subroutine linear_transform_backward_nr(qprop, xprop,bead)
      implicit none
      double precision::        xprop(:), qprop(:)
      double precision, allocatable:: xprop_nr(:)
      integer::                 i,j,bead
      allocate(xprop_nr(1:n+1))
      if (bead .gt. 0) then
         xprop_nr(2:n+1)=qprop(1:n) + beadvec(1:n, bead)
      else 
         xprop_nr(2:n+1)=qprop(1:n)
      end if
      ! call sinft(xprop_nr)
      xprop(1:n)= xprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
      deallocate(xprop_nr)
      return
    end subroutine linear_transform_backward_nr

    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode forward transformation for linear polymer
    subroutine ring_transform_forward_nr(xprop, qprop)
      implicit none
      double precision, intent(in):: xprop(:)
      double precision, intent(out):: qprop(:)
      double precision:: qprop_nr(n)
      integer::                 i,j,k, ier

      qprop_nr(:)= xprop(:)

      ier = DftiComputeForward(fft_handle, qprop_nr)

      qprop(1)= qprop_nr(1)
      qprop(n)= qprop_nr(n)
      do i=2,n-1,2
         qprop(i)= qprop_nr(i)*sqrt(2.0d0)
         qprop(i+1)= -qprop_nr(i+1)*sqrt(2.0d0)
      end do
      return
    end subroutine ring_transform_forward_nr
    !-----------------------------------------------------
    !-----------------------------------------------------
    !normal mode backward transformation for linear polymer
    subroutine ring_transform_backward_nr(qprop, xprop)
      implicit none
      double precision, intent(in):: qprop(:)
      double precision, intent(out):: xprop(:)
      double precision:: xprop_nr(n)
      integer::                 i,j,k,ier

      xprop_nr(1)= qprop(1)
      xprop_nr(n)= qprop(n)
      do i=2,n-1,2
         xprop_nr(i)=qprop(i)/sqrt(2.0d0)
         xprop_nr(i+1)=-qprop(i+1)/sqrt(2.0d0)
      end do
      ier = DftiComputeBackward(fft_handle, xprop_nr)
      xprop(:)= xprop_nr(:)
      return
    end subroutine ring_transform_backward_nr

    !TODO: replace idof calculations with this function
    function calcidof(i,j) !i for dimensions, j for atom
      implicit none
      integer:: calcidof
      integer, intent(in):: i,j
      calcidof= (j-1)*ndim + i
    end function calcidof

    function centroid(x)
      implicit none
      double precision, intent(in):: x(:)
      double precision:: centroid
      integer:: i

      do i=1,n
         centroid=centroid+ x(i)/dble(n)
      end do
      return
    end function centroid


    subroutine harmonicsampling(x,p, potvals)
      implicit none
      double precision, intent(inout):: x(:,:,:), p(:,:,:)
      double precision, intent(out):: potvals
      double precision,allocatable:: rp(:), tempx(:), pk(:), rk(:,:)
      double precision:: lambda, stdev
      integer:: i,j,k, idof

      allocate(rp(n),tempx(totdof), pk(n))
      allocate(rk(n,ndof))
      stdev= 1.0d0/sqrt(betan)
      x(:,:,:)= 0.0d0
      p(:,:,:)= 0.0d0
      potvals=0.0d0

      do i=1,ndim
         do j=1,natom
            idof= calcidof(i,j)
            if (n .gt. 1) then
               !---------------------------
               !generate random numbers for chain
               errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,rp,0.0d0,stdev)!ringpolymer
               !---------------------------
               !scale them to have the right standard deviation
               do k=1, n
                  lambda=sqrt(transfreqs(idof) + lam(k)**2)
                  rp(k)= rp(k)/(lambda*sqrt(mass(j)))
                  potvals=potvals + 0.5d0*mass(j)*lambda**2*rp(k)**2
               end do
               !---------------------------
               !transform to site-local cartesian
               if (ring) then
                  if (use_fft) then
                     call ring_transform_backward_nr(rp, rk(:,idof))
                  else
                     call ring_transform_backward(rp, rk(:,idof))
                  end if
               else
                  call linear_transform_backward(rp, rk(:,idof), idof)
               end if
            else
               errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,1,rk(1,:),0.0d0,stdev)!ringpolymer
               rk(1,idof)= rk(1,idof)/(sqrt(mass(j)*transfreqs(j)))
            end if
            errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,p(:,i,j),0.0d0,stdev)!momenta
            p(:,i,j)= p(:,i,j)*sqrt(mass(j))
         end do
      end do
      !---------------------------
      !transform to global cartesian coordinates
      do k=1,n
         x(k,:,:)= reshape(matmul(transpose(hess),rk(k,:)), (/ndim,natom/))
         x(k,:,:)= x(k,:,:) + transition(:,:)
      end do

      deallocate(tempx, rp, rk, pk)
      
      return
    end subroutine harmonicsampling

    
  end module general
