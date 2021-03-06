#TODO: run checks for correct folders created in makefile and stop tracking the folders altogether
###################################################################################
#Define all the compilers and their options
FC= ifort
MPIFC= mpiifort
FFLAGS= -warn -g -O2 -i8 -r8 -w -no-wrap-margin -module ../modules
FLIBS=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
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

#1D eckart barrier :
ECKARTFILES= $(OBJDIR)/potential_eckart.o

#1D harmonic lattice:
HARMFILES= $(OBJDIR)/potential_harmoniclattice.o

#1D harmonic potential:
SHOFILES= $(OBJDIR)/potential_SHO.o

#1D quartic potential:
QUARTFILES= $(OBJDIR)/potential_quartic.o

rpmd_eckart_rate: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_rates.o $(OBJDIR)/verletmodule.o

rpmd_eckart_flux: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_fluxflux.o $(OBJDIR)/verletmodule.o

rpmd_harmonic_heat_ser: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_heat.o $(OBJDIR)/verletmodule.o

rpmd_harmonic_heat_par: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_heat.o $(OBJDIR)/verletmodule.o

rpmd_sho_position_ser: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_position.o $(OBJDIR)/verletmodule.o

rpmd_sho_position_par: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_position.o $(OBJDIR)/verletmodule.o

rpmd_quartic_position_ser: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_position.o $(OBJDIR)/verletmodule.o

rpmd_quartic_position_par: POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/estimators_position.o $(OBJDIR)/verletmodule.o

###################################################################################
#Compilation commands for object files
$(OBJDIR)/%.o: %.f
	$(FC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/heattransport_par.o: heattransport_par.f90
	$(MPIFC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/position_par.o: position_par.f90
	$(MPIFC) $(INCLUDE) -c $(FFLAGS) $< -o $@

$(OBJDIR)/mkl_dfti.o: ${MKLROOT}/include/mkl_dfti.f90
	$(FC) -c $(FFLAGS) $< -o $@

###################################################################################
#Rules for the final executables

################################
#eckart barrier potential:
rpmd_eckart_rate: $(PREPOTFILES) $(ECKARTFILES) $$(POSTPOTFILES) $(OBJDIR)/rpmd.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ECKARTFILES) $(POSTPOTFILES) $(OBJDIR)/rpmd.o -o $(BUILDDIR)/$@ $(FLIBS)

rpmd_eckart_flux: $(PREPOTFILES) $(ECKARTFILES) $$(POSTPOTFILES) $(OBJDIR)/fluxflux.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ECKARTFILES) $(POSTPOTFILES) $(OBJDIR)/fluxflux.o -o $(BUILDDIR)/$@ $(FLIBS)


################################
#harmonic lattice potential:
rpmd_harmonic_heat_ser: $(PREPOTFILES) $(HARMFILES) $$(POSTPOTFILES) $(OBJDIR)/heattransport_ser.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(HARMFILES) $(POSTPOTFILES) $(OBJDIR)/heattransport_ser.o -o $(BUILDDIR)/$@ $(FLIBS)

rpmd_harmonic_heat_par: $(PREPOTFILES) $(HARMFILES) $$(POSTPOTFILES) $(OBJDIR)/heattransport_par.o
	$(MPIFC) $(FFLAGS) $(PREPOTFILES) $(HARMFILES) $(POSTPOTFILES) $(OBJDIR)/heattransport_par.o -o $(BUILDDIR)/$@ $(FLIBS)

################################
#harmonic potential:
rpmd_sho_position_ser: $(PREPOTFILES) $(SHOFILES) $$(POSTPOTFILES) $(OBJDIR)/position_ser.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(SHOFILES) $(POSTPOTFILES) $(OBJDIR)/position_ser.o -o $(BUILDDIR)/$@ $(FLIBS)

rpmd_sho_position_par: $(PREPOTFILES) $(SHOFILES) $$(POSTPOTFILES) $(OBJDIR)/position_par.o
	$(MPIFC) $(FFLAGS) $(PREPOTFILES) $(SHOFILES) $(POSTPOTFILES) $(OBJDIR)/position_par.o -o $(BUILDDIR)/$@ $(FLIBS)

################################
#quartic potential:
rpmd_quartic_position_ser: $(PREPOTFILES) $(QUARTFILES) $$(POSTPOTFILES) $(OBJDIR)/position_ser.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(QUARTFILES) $(POSTPOTFILES) $(OBJDIR)/position_ser.o -o $(BUILDDIR)/$@ $(FLIBS)

rpmd_quartic_position_par: $(PREPOTFILES) $(QUARTFILES) $$(POSTPOTFILES) $(OBJDIR)/position_par.o
	$(MPIFC) $(FFLAGS) $(PREPOTFILES) $(QUARTFILES) $(POSTPOTFILES) $(OBJDIR)/position_par.o -o $(BUILDDIR)/$@ $(FLIBS)

###################################################################################
#Rules for cleanup

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(MODDIR)/*.f90
