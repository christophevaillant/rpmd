###################################################################################
#Define all the compilers and their options
FC= ifort
MPIFC= mpif90
FFLAGS= -warn -g -O2 -i8 -r8 -w -no-wrap-margin -module ../modules
FLIBS_PAR=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
FLIBS_SEQ=-L${MKLROOT}/lib/intel64 -mkl=sequential -lmkl_rt -lpthread -lm -ldl
INCLUDE= -I${MKLROOT}/include
BUILDDIR= ../bin
OBJDIR= ../objects
MODDIR= ../modules
###################################################################################
#List the object files and their dependencies
#Define some variables with common files
PREPOTFILES= $(OBJDIR)/nr_fft.o $(OBJDIR)/general.o
POSTPOTFILES= $(OBJDIR)/instantonmod.o $(OBJDIR)/verletmodule.o
RATEFILES= $(OBJDIR)/estimators_rates.o

#1D eckart barrier :
ECKARTFILES= $(OBJDIR)/potential_eckart.o

###################################################################################
#Compilation commands for object files
$(OBJDIR)/%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@


###################################################################################
#Rules for the final executables

################################
#eckart barrier potential:
rpmd_eckart_rate: $(PREPOTFILES) $(ECKARTFILES) $(RATEFILES) $(POSTPOTFILES) $(OBJDIR)/rpmd.o
	$(FC) $(FFLAGS) $(PREPOTFILES) $(ECKARTFILES) $(RATEFILES) $(POSTPOTFILES) $(OBJDIR)/rpmd.o -o $(BUILDDIR)/$@ $(FLIBS_SEQ)


###################################################################################
#Rules for cleanup

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(MODDIR)/*.f90
