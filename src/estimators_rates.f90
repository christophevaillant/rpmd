module estimators
implicit none

integer::  nestim, ntime
logical, allocatable:: whichestim(:)
double precision, allocatable:: tcf(:,:)

contains

subroutine init_estimators()
  implicit none

  ntime= mod(NMC,noutput)
  allocate(whichestim(nestim), tcf(nestim,ntime))

end subroutine init_estimators


end module estimators
