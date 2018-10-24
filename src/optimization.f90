module optimization
use potential_derivatives
use variables
use constants
use classes
use system_general
implicit none

public:: instanton, ringpolymer_pot, ringpolymer_hess, gtol

private

integer:: totdof, ndof
real(adequate):: betan

!TODO: implement Hessian update approximation
  
contains
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine streambed(x, force, modes, freqs, order, nzero, maxstep)
    !Stream-bed walking algorithm (Nicholls JCP 92, pp 340, 1990), as
    !adapted by J O Richardson for instantons.

    !Inputs:
    !x(ndof) (real): contains array of variables to be optimized
    !force(ndof) (real): vector of derivatives of potential
    !modes(ndof,ndof) (real): eigenvectors of second derivative matrix of potential
    !freqs(ndof) (real): eigenvalues of second derivative matrix of potential
    !order (int): =0 for minimization, =1 for max, =2 for TS search
    !nzero (int): number of hessian zeros (i.e. conserved motions) to skip
    !maxstep (real): largest step

    implicit none
    integer, intent(in)::  order, nzero
    real(adequate), intent(inout):: x(:)
    real(adequate), intent(in):: maxstep
    real(adequate), intent(in):: force(:), modes(:,:), freqs(:)
    real(adequate), allocatable::  steps(:), fake(:), normalforce(:), cartsteps(:)
    real(adequate):: energy, lam, alpha, normstep
    integer:: i, maxforceindex

    allocate(fake(totdof), cartsteps(totdof))
    allocate(normalforce(totdof), steps(totdof))

    !convert force to normal mode coords
    normalforce(:)=0.0
    normalforce(:)= matmul(force,modes)
    steps(:)=0.0d0

    !take steps according to different criteria (difference between min and TS)
    !----------------------------------------
    if (order.lt.2) then
       !perform minimization/maximization
       if (freqs(1) .lt. 0.0d0) then
          lam= freqs(1) - abs(normalforce(1)/maxstep)
          steps(:)= normalforce(:)/(lam-freqs(:))
       else
          steps(:)= normalforce(:)/(freqs(:))
       end if
       if (order.eq.1) steps(:)=-steps(:)

    !----------------------------------------
    else if (order .eq. 2) then
       !perform transition state search
       if (freqs(1).gt. 0.0) then 
          !all frequencies positive because they're sorted
          if (0.5*freqs(2+nzero) .gt. freqs(1)) then
             alpha=1.0
             lam= (2.0*freqs(1) + freqs(2+nzero))/4.0
          else
             alpha= (freqs(2+nzero) - freqs(1))/freqs(2+nzero)
             lam= 0.25d0*(3.0d0*freqs(1)+ freqs(2+nzero))
          end if
       else if (freqs(2+nzero) .lt. 0.0) then
          !first 2 freqs negative
          !here lies JOR's modifications for ring polymer stuff
          if (freqs(2+nzero) .ge. freqs(1)/2.0) then
             alpha=1.0d0
             lam= 0.25d0*(freqs(1) + 2.0d0*freqs(2+nzero))
          else
             alpha= (freqs(1) - freqs(2+nzero))/freqs(2+nzero)
             lam= 0.25d0*(freqs(1) + 3.0d0*freqs(2+nzero))
          end if
       else 
          !first freq negative, rest are positive
          alpha=1.0
          lam= 0.25d0*(freqs(1) + freqs(2+nzero))
       end if
       steps(:)= alpha*normalforce(:)/(lam-freqs(:))
    end if

    !rescale the step to be below max
    normstep= norm2(steps)
    if (normstep .gt. maxstep) steps(:)=steps(:)*maxstep/normstep
    !convert back to Cartesians
    cartsteps= matmul(modes, steps)
    do i=1, totdof
       x(i)= x(i)+ cartsteps(i)
    end do
 

    deallocate(steps, normalforce, fake, cartsteps)
  end subroutine streambed
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!Routine to approximately update the Hessian

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!Driver Routine to calculate a semiclassical instanton configuration
  subroutine instanton(polymer, maxstep, maxiter, constrain, a, b, approxhess, gtol)
    implicit none
    type(T_ring_polymer), intent(inout):: polymer
    integer, intent(in):: maxiter
    real(adequate), intent(in):: maxstep, a(:, :), b(:, :), gtol
    logical, intent(in):: constrain, approxhess
    real(adequate), allocatable:: coords(:), force(:), hess(:,:)
    real(adequate), allocatable:: modes(:,:), freqs(:),work(:)
    real(adequate):: gtol, rms
    integer, allocatable::  iwork(:)
    integer:: i,j,k, idof, info, lwork, liwork, ii, minindex
    
    !-------------------------
    !Inputs:
    !polymer (T_ring_polymer): the initial guess
    !maxstep (real): the largest allowed step during optimization
    !maxiter (int): maximum number of iterations
    !constrain (logical): whether or not to fix end points
    !a/b(NDimen, NAtoms) (real): end points if constrain is true
    !approxhess (logical): whether or not to approximate the hessian (NOT IMPLEMENTED)
    !gtol (real): sets the tolerance for the gradient
    
    betan=beta/real(NBeads)
    write(*,*) "Betan=", betan, beta
    totdof= Nbeads*Natoms*NDimen
    ndof= Natoms*NDimen

    !allocate various arrays
    lwork=1+ 6*totdof + 2*totdof**2
    liwork= 3 + 5*totdof
    allocate(coords(totdof),force(totdof),hess(totdof, totdof))
    allocate(modes(totdof,totdof),freqs(totdof))
    allocate(work(lwork), iwork(liwork))
    hess(:,:)=0.0
    force(:)=0.0
    freqs(:)=0.0
    modes(:,:)=0.0
    do ii=1, maxiter
       hess(:,:)=0.0
       force(:)=0.0
       freqs(:)=0.0
       modes(:,:)=0.0
       call ringpolymer_grad(polymer, constrain,a,b,force)

       !calculate RMS & test finishing criteria
       rms= sqrt(dot_product(force,force)/real(totdof))
       if (RMS .lt. gtol) exit

       !calculate normal mode transforms
       call ringpolymer_hess(polymer, hess)
       call dsyevd('V', 'U', totdof, hess, totdof, freqs,work,lwork,iwork,&
            liwork,info)
       modes(:,:)= hess(:,:)
       minindex= minloc(abs(freqs),dim=1)
       if (modes(minindex,minindex) .lt. 0.0) modes(:,:)= -modes(:,:)

       !copy out the ring polymer coordinates into a flat array
       do i=1, NBeads
          do j=1, Ndimen
             do k= 1,Natoms
                idof= i + NBeads*((j-1) + Ndimen*(k-1))
                coords(idof)= polymer%beads(i)%coords(j,k)
             end do
          end do
       end do

       !optimize step
       if (constrain) then
          call streambed(coords, force, modes, freqs, 2, gtol, 1, maxstep)
       else
          call streambed(coords, force, modes, freqs, 2, gtol, 0, maxstep)
       end if

       !convert back to polymer
       do i=1, NBeads
          do j=1, Ndimen
             do k= 1,Natoms
                idof= i + NBeads*((j-1) + Ndimen*(k-1))
                polymer%beads(i)%coords(j,k)= coords(idof)
             end do
          end do
       end do

    end do
    if (ii.ge.maxiter+1) then
         write(*,*) "WARNING: Maximum number of iterations reached with no convergence."
         write(*,*) "WARNING: final RMS was:", RMS
      end if
    deallocate(coords, force, hess, modes, freqs, work, iwork)

  end subroutine instanton

  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !Function to calculate energy from ring polymer
  function ringpolymer_pot(polymer,constrain,a,b)
    implicit none
    type(T_ring_polymer), intent(inout):: polymer
    logical, intent(in)::  constrain
    real(adequate), intent(in)::  a(:,:), b(:,:)
    integer::            i,j,k
    double precision::   ringpolymer_pot, UM

    UM=0.0d0
    do i=1, Nbeads-1, 1
       UM=UM+ energy(polymer%beads(i)%coords)
       do j=1, NDimen
          do k=1, NAtoms
             UM=UM+ (0.5d0*mass(k)/betan**2)*(polymer%beads(i+1)%coords(j,k)-&
                  polymer%beads(i)%coords(j,k))**2
          end do
       end do
    end do
    if (constrain) then
       UM=UM+ energy(a(:,:))+ energy(b(:,:))+ energy(polymer%beads(Nbeads)%coords)
       do j=1, NDimen
          do k=1, NAtoms
             UM=UM+ (0.5d0*mass(k)/betan**2)*(polymer%beads(1)%coords(j,k)-a(j,k))**2
             UM=UM+ (0.5d0*mass(k)/betan**2)*(b(j,k)-polymer%beads(NBeads)%coords(j,k))**2
          end do
       end do
    else
       UM=UM+ energy(polymer%beads(Nbeads)%coords)
       do j=1, NDimen
          do k=1, NAtoms
             UM=UM+ (0.5d0*mass(k)/betan**2)*(polymer%beads(1)%coords(j,k)-&
                  polymer%beads(NBeads)%coords(j,k))**2
          end do
       end do
    end if

    ringpolymer_pot= UM

    return
  end function ringpolymer_pot
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!Function to calculate force from ring polymer
subroutine ringpolymer_grad(polymer, constrain, a, b, grad)
  implicit none
  type(T_ring_polymer), intent(inout):: polymer
  logical, intent(in)::  constrain
  real(adequate), intent(in)::  a(:,:), b(:,:)
  real(adequate), intent(out):: grad(:)
  real(adequate), allocatable:: tempcoords(:)
  integer:: i,j,k,idof

  allocate(tempcoords(totdof))
  grad(:)=0.0d0
    do i=1, Nbeads
       call get_en_force(polymer%beads(i)%coords,polymer%beads(i)%pot_val,&
            polymer%beads(i)%force_vals)
       do j=1,NDimen
          do k=1,NAtoms
             idof= i + NBeads*((j-1) + Ndimen*(k-1))
             if (i.eq.1) then
                if (constrain) then
                   grad(idof)=mass(k)*(2.0*polymer%beads(1)%coords(j,k)&
                        - a(j,k) - polymer%beads(2)%coords(j,k))/betan**2
                else
                   grad(idof)=mass(k)*(2.0*polymer%beads(1)%coords(j,k) - &
                        polymer%beads(2)%coords(j,k)- polymer%beads(NBeads)%coords(j,k))/betan**2
                end if

             else if (i.eq.NBeads) then
                if (constrain) then
                   grad(idof)=(mass(k)*2.0*polymer%beads(NBeads)%coords(j,k)&
                        - polymer%beads(NBeads-1)%coords(j,k) - b(j,k))/betan**2
                else
                   grad(idof)=mass(k)*(2.0*polymer%beads(NBeads)%coords(j,k) &
                        - polymer%beads(NBeads-1)%coords(j,k)- polymer%beads(1)%coords(j,k))/betan**2
                end if
             else
                grad(idof)= mass(k)*(2.0*polymer%beads(i)%coords(j,k)&
                     - polymer%beads(i-1)%coords(j,k) - polymer%beads(i+1)%coords(j,k))/betan**2
             end if

             grad(idof)= grad(idof)- polymer%beads(i)%force_vals(j,k)
             grad(idof)= grad(idof)/mass(k)
          end do
       end do
    end do

end subroutine ringpolymer_grad

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!Function to calculate Hessian from linear polymer
subroutine ringpolymer_hess(polymer, hess)
  implicit none
  type(T_ring_polymer), intent(in):: polymer
  real(adequate), intent(out):: hess(:,:)
  real(adequate), allocatable:: hessbead(:,:), dummy(:)
  integer:: i1,i2,j1,k1,idof1,j2,k2,idof2
  integer:: fulldof1, fulldof2, index

  allocate(hessbead(ndof,ndof), dummy(ndof))
    
  hess=0.0d0
  do i1=1, NBeads, 1
     hessbead=0.0d0
     call hessian_find(polymer%beads(i1)%coords, hessbead, dummy, ndof)
     do j1=1,NDimen
        do k1=1,NAtoms
           do i2=1, Nbeads
           do j2=1,NDimen
              do k2=1,NAtoms
                 !The DoF label for the individual bead
                 idof1= (k1-1)*NDimen + j1
                 idof2= (k2-1)*NDimen + j2
                 !The DoF label for the whole matrix
                 fulldof1=ndof*(i1-1) + idof1
                 fulldof2=ndof*(i2-1) + idof2

                 if (i1.eq.i2) then
                    hess(fulldof1,fulldof2)=hessbead(idof1, idof2)/sqrt(mass(k1)*mass(k2))
                    if (idof1.eq.idof2) hess(fulldof1,fulldof2)= &
                         hess(fulldof1,fulldof2)+2.0d0/betan**2
                 else if ((i1 .eq. i2+1) .or. (i1 .eq. i2-1)) then
                    if (idof1.eq.idof2) hess(fulldof1,fulldof2)=-1.0/betan**2
                 end if
                 if ((i1.eq.NBeads) .and. (i2.eq.1)) then
                    if (idof1 .eq. idof2) then
                    hess(fulldof1,fulldof2)=-1.0/betan**2
                    hess(fulldof2,fulldof1)=-1.0/betan**2
                 end if
                 end if
              end do
           end do
        end do
     end do
  end do
  end do
  deallocate(hessbead, dummy)
  return

end subroutine ringpolymer_hess

subroutine polymer_to_flat(polymer, coords)
implicit none
integer:: i,j,k,idof
real(adequate), intent(out):: coords(:)
type(T_ring_polymer), intent(in):: polymer

       do i=1, NBeads
          do j=1, Ndimen
             do k= 1,Natoms
                idof= i + NBeads*((j-1) + Ndimen*(k-1))
                coords(idof)= polymer%beads(i)%coords(j,k)
             end do
          end do
       end do

end subroutine polymer_to_flat

end module optimization
