module estimators
  use general

implicit none

logical, allocatable:: whichestim(:)
double precision, allocatable:: transition(:,:), normalvec(:,:), transfreqs(:)
double precision, allocatable:: lambda(:,:), hess(:,:)

contains

subroutine init_estimators()
  implicit none
  integer:: i,j,k, idof
  !TODO: estimatormod- implement list of estimators

  ntime= mod(NMC,noutput)
  allocate(whichestim(nestim), tcf(nestim,ntime))
  allocate(lambda(ndof,n))
  
  do i=1, ndim
     do j=1, natom
        do k=1, n
           idof= calcidof(i,j)
           lambda(idof, k)= 4.0*((omegan**2 - transfreqs(idof)**2) &
                - omegan**2*cos(2.0*pi*k/n))
        end do
     end do
  end do

end subroutine init_estimators

subroutine finalize_estimators()
  implicit none
  deallocate(whichestim,tcf)
end subroutine finalize_estimators

subroutine estimator(x,p,tcf)
  implicit none
  double precision, intent(in):: x(:,:,:), p(:,:,:)
  double precision, intent(out):: tcf
  double precision,allocatable:: v(:,:)
  double precision:: s
  integer:: i,j

  allocate(v(ndim,natom))
  do i=1,ndim
     do j=1,natom
        v(i,j)= centroid(x(:,i,j)) - transition(i,j)
     end do
  end do
  s= dot_product(reshape(v,(/ndof/)), reshape(normalvec,(/ndof/)))
  if (s .gt. 0) then
     tcf=1.0d0
  else
     tcf=0.0d0
  end if

  deallocate(v)
  return
end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, weight)
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision, intent(out)::  weight
    double precision, allocatable::   vel(:), tempp(:), pos(:), tempx(:)
    double precision, allocatable::  rk(:,:)
    double precision::                stdev
    integer::                         i,j,k, idof, imin(1)

    allocate(tempp(n), vel(n),pos(n), tempx(n))
    allocate(rk(n,ndof))
    do i=1, ndim
       do j=1, natom
          !---------------------------
          !generate random numbers
          idof= calcidof(i,j)
          stdev= 1.0d0/sqrt(betan)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,pos,0.0d0,stdev)
          errcode_normal = vdrnggaussian(rmethod_normal,stream_normal,n,vel,0.0d0,stdev)
          !---------------------------
          !scale them to have the right standard deviation
          do k=1, n
             if (k.eq.1) then
                pos(1)=0.0d0 !set centroid mode to 0.0
             else
                pos(k)= pos(k)/sqrt(mass(k)*lambda(idof,k))
             end if
             vel(k)= vel(k)*sqrt(beadmass(j,k))
          end do
          !---------------------------
          !transform to non-ring polymer normal mode coords
          if (ring) then
             call ring_transform_backward(pos, tempx)
             call ring_transform_backward(vel, tempp)
          else
             call linear_transform_backward(pos, tempx, idof)
             call linear_transform_backward(vel, tempp, 0)
          end if
          rk(:,idof)= tempx(:)
          p(:,i,j)= tempp(:)
       end do
    end do
    !---------------------------
    !transform to cartesian coordinates
    do k=1,n
       x(k,:,:)= reshape(matmul(transpose(hess),rk(k,:)),(/ndim,natom/))&
            - transition(:,:)
    end do
    deallocate(vel, tempp,pos, tempx)
    return
  end subroutine init_path

end module estimators
