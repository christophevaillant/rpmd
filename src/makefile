#TODO: run checks for correct folders created in makefile and stop tracking the folders altogether
###################################################################################
#Define all the compilers and their options
FC= ifort
MPIFC= mpif90
FFLAGS= -warn -g -O2 -i8 -r8 -w -no-wrap-margin -module ../modules
FLIBS_PAR=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
FLIBS_SEQ=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl -L. -lfftpack -lgfortran
INCLUDE= -I${MKLROOT}/include
BUILDDIR= ../bin
OBJDIR= ../objects
MODDIR= ../modules
###################################################################################
#List the object files and their dependencies
#Define some variables with common files
PREPOTFILES= $(OBJDIR)/mkl_dfti.o $(OBJDIR)/general.o
.SECONDEXPANSION:
POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/verletmodule.o

# RATEFILES= $(OBJDIR)/estimators_rates.o

#1D eckart barrier :
ECKARTFILES= $(OBJDIR)/potential_eckart.o

#1D asymmetric barrier :
ASYMMFILES= $(OBJDIR)/potential_asymmetric.o

rpmd_eckart_rate: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_rates.o $(OBJDIR)/verletmodule.o

rpmd_eckart_flux: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_fluxflux.o $(OBJDIR)/verletmodule.o

rpmd_asymm_rate: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_rates.o $(OBJDIR)/verletmodule.o

rpmd_asymm_flux: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_fluxflux.o $(OBJDIR)/verletmodule.o


###################################################################################
#Compilation commands for object files
$(OBJDIR)/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/mkl_dfti.o: ${MKLROOT}/include/mkl_dfti.f90
	$(FC) -c $(FFLAGS) $< -o $@

###################################################################################
#Rules for the final executables

################################
#eckart barrier potential:
rpmd_eckart_rate: $(PREPOTFILES) $(ECKARTFILES) $$(POSTPOTFILES) $(OBJDIR)/rpmd.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ECKARTFILES) $(POSTPOTFILES) $(OBJDIR)/rpmd.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpmd_eckart_flux: $(PREPOTFILES) $(ECKARTFILES) $$(POSTPOTFILES) $(OBJDIR)/fluxflux.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ECKARTFILES) $(POSTPOTFILES) $(OBJDIR)/fluxflux.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

################################
#asymmetric eckart barrier potential:
rpmd_asymm_rate: $(PREPOTFILES) $(ASYMMFILES) $$(POSTPOTFILES) $(OBJDIR)/rpmd.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ASYMMFILES) $(POSTPOTFILES) $(OBJDIR)/rpmd.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

rpmd_asymm_flux: $(PREPOTFILES) $(ASYMMFILES) $$(POSTPOTFILES) $(OBJDIR)/fluxflux.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ASYMMFILES) $(POSTPOTFILES) $(OBJDIR)/fluxflux.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)

###################################################################################
#Rules for cleanup

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(MODDIR)/*.f90
