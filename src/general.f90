include 'mkl_vsl.fi'
module general
  use MKL_VSL_TYPE
  use MKL_VSL
  use nr_fft

  implicit none

  integer::                         NMC, Noutput,nestim, ntime
  integer::                         errcode_poisson, rmethod_poisson
  integer::                         brng_poisson, seed_poisson
  type (vsl_stream_state)::         stream_poisson,stream_normal
  integer::                         errcode_normal, rmethod_normal, brng_normal
  integer::                         seed_normal
  integer::                         n, ndim, ndof, natom, xunit, totdof
  double precision, allocatable::   tcf(:,:)
  double precision,allocatable::    transmatrix(:,:),beadvec(:,:)
  double precision,allocatable::    beadmass(:,:), lam(:)
  logical::                         ring
  double precision, parameter::    pi=3.14159265358979d0
  double precision::               beta, betan, UMtilde, omegan
  double precision, allocatable::  well1(:,:), well2(:,:), mass(:)
  character, allocatable::         label(:)
  logical::                        fixedends

contains

  !-----------------------------------------------------
  !-----------------------------------------------------
  !allocate normal mode arrays

  subroutine alloc_nm(iproc)
    implicit none
    integer::    itime, irate, imax
    integer, intent(in):: iproc

    allocate(transmatrix(n,n), beadmass(natom,n), lam(n))
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
    deallocate(transmatrix,beadmass,lam)
    if (.not. ring) deallocate(beadvec)
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
    integer::                 i,j

    call dsymv('T', n, 1.0d0,transmatrix, n, xprop,1,0.0d0, qprop,1)

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

    call dsymv('U', n, 1.0d0, transmatrix, n, qprop,1,0.0d0, xprop,1)

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
    call sinft(qprop_nr)
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
    call sinft(xprop_nr)
    xprop(1:n)= xprop_nr(2:n+1)*sqrt(2.0d0/(dble(n+1)))
    deallocate(xprop_nr)
    return
  end subroutine linear_transform_backward_nr

  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode forward transformation for linear polymer
  subroutine ring_transform_forward_nr(xprop, qprop)
    implicit none
    double precision::        xprop(:), qprop(:)
    double precision, allocatable:: qprop_nr(:)
    integer::                 i,j
    !TODO: implement numerical recipes FFT

    write(*,*) "Not implemented!"
    stop
    return
  end subroutine ring_transform_forward_nr
  !-----------------------------------------------------
  !-----------------------------------------------------
  !normal mode backward transformation for linear polymer
  subroutine ring_transform_backward_nr(qprop, xprop)
    implicit none
    double precision::        xprop(:), qprop(:)
    double precision, allocatable:: xprop_nr(:)
    integer::                 i,j
    !TODO: implement numerical recipes inverse FFT

    write(*,*) "Not implemented!"
    stop
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
end module general
