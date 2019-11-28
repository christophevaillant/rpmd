module estimators
  use general
  use instantonmod
  use potential

  implicit none

  public:: init_estimators, estimator, finalize_estimators,init_path, whichestim, convection

  private

  logical, allocatable:: whichestim(:)
  logical:: convection

contains

  subroutine init_estimators()
    implicit none
    integer:: i,j,k,kk, idof
    !TODO: estimatormod- implement list of estimators

    if (outputfbar) open(200,file="fbar.dat")

    ntime= NMC/Noutput
    allocate(whichestim(nestim))

  end subroutine init_estimators

  subroutine finalize_estimators()
    implicit none
    deallocate(whichestim)
  end subroutine finalize_estimators

  subroutine estimator(x,p,tcfval)
    implicit none
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision, intent(out):: tcfval(:)
    double precision:: s, energy, prob, stdev, q1, q2
    double precision:: centroidxi, centroidpi, centroidxj, centroidpj
    integer:: i,j,k, idof

    tcfval(1)=0.0d0
    do i=1, natom
       centroidxi= centroid(x(:,1,i))
       q1=centroidxi - lattice(1,i)
       centroidpi= centroid(p(:,1,i))
       if (i .lt. natom) then
          centroidxj= centroid(x(:,1,i+1))
          q2= centroidxj - lattice(1,i+1)
          centroidpj= centroid(p(:,1,i+1))
          tcfval(1)= tcfval(1) + &
               0.5d0*(centroidpi + centroidpj)*interharm**2*(centroidxj-centroidxi)*(q1-q2)
       end if
       if (i .gt. 1) then
          centroidxj= centroid(x(:,1,i-1))
          q2= centroidxj - lattice(1,i-1)
          centroidpj= centroid(p(:,1,i-1))
          tcfval(1)= tcfval(1) + &
               0.5d0*(centroidpi + centroidpj)*interharm**2*(centroidxi-centroidxj)*(q2-q1)
       end if
    end do

    ! tcfval(1)=0.0d0
    ! do i=1,natom
    !    energy=0.0d0
    !    do j=1,n
    !       if (i .lt. natom) energy= energy + &
    !            0.5d0*(x(j,1,i+1) - x(j,1,i))*interforce(x,i,j)*(p(j,1,i) + p(j,1,i+1))/mass(i)
    !       if (i .gt. 1) energy= energy + &
    !            0.5d0*(x(j,1,i) - x(j,1,i-1))*interforce(x,i-1,j)*(p(j,1,i-1) + p(j,1,i))/mass(i)
    !       if (convection) energy= energy+ p(j,1,i)*siteenergy(p,x,i,j)/mass(i)
    !    end do
    !    tcfval(1)= tcfval(1)+energy/dble(N)
    ! end do

    tcfval(2)= tcfval(1)

    return
  end subroutine estimator

  !-----------------------------------------------------
  !-----------------------------------------------------
  !function for initializing a path
  subroutine init_path(x, p, factors,weight)
    double precision, intent(inout):: x(:,:,:), p(:,:,:)
    double precision, intent(out)::   weight,factors(:)
    double precision::              potvals, ringpot, potdiff

    call harmonicsampling(x, p, potvals)
    ringpot= UM(x)
    potdiff= (ringpot- potvals)
    weight= exp(-betan*potdiff)

    !work out initial current
    call estimator(x,p,factors)
    factors(2)=1.0d0
    
    return
  end subroutine init_path


end module estimators
