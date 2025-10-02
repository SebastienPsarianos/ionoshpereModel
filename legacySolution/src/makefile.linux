#
# Makefile for mglin & mgfas: Linear and Full Approximation Storage (FAS)
#                             Full Multigrid Solvers.
# also
# Makefile for 2D Ionospheric Electric Potential Model.
#
# Linux Version
#
SOURCEDIR = ~/groth/Ionosphere/src
OBJ1 = mglin.o
OBJ2 = mgfas.o
OBJ3 = ionosphere.o
IONOBJ = IONOSPHERE_conductance.o IONOSPHERE_multigrid.o \
         IONOSPHERE_fgrid.o IONOSPHERE_cgrid.o \
         IONOSPHERE_restrict.o IONOSPHERE_prolong.o \
         IONOSPHERE_smoother.o IONOSPHERE_residual.o \
         IONOSPHERE_update.o IONOSPHERE_matcopy.o \
         IONOSPHERE_matzero.o IONOSPHERE_util.o \
         IONOSPHERE_magdim.o IONOSPHERE_magfac.o \
         IONOSPHERE_magbfield.o IONOSPHERE_magvel.o IONOSPHERE_current.o
CPTGOBJ = CPTG_util.o
NROBJ = NR_mglin.o NR_mgfas.o \
        NR_addint.o NR_interp.o NR_lop.o NR_relax.o NR_relax2.o \
        NR_resid.o NR_rstrct.o NR_slvsml.o NR_slvsm2.o \
        NR_copy.o NR_matadd.o NR_matsub.o NR_fill0.o NR_anorm2.o \
        NR_util.o
IONINC = IONOSPHERE.h IONOSPHERE_util.h
CPTGINC = CPTG_util.h
NRINC = NR.h NR_util.h
LIB = -lm
EXE1 = mglin.exe
EXE2 = mgfas.exe
EXE3 = ionosphere.exe
CC = gcc
BIND = gcc
CPU = 
OPT = -O
DEBUG =
#
# Define the dependencies for the three main programs.
#
all : $(EXE1) $(EXE2) $(EXE3)
	@echo Programs MGLIN, MGFAS, and IONOSPHERE have been brought up to date.

mglin : $(EXE1)
	@echo Program MGLIN has been brought up to date.

mgfas : $(EXE2)
	@echo Program MGFAS has been brought up to date.

ionosphere : $(EXE3)
	@echo Program IONOSPHERE has been brought up to date.
#
# Define the dependencies for the three main executables.
#
$(EXE1) : $(OBJ1) $(NROBJ)
	$(BIND) -o $(EXE1) $(OBJ1) $(NROBJ) $(LIB)

$(EXE2) : $(OBJ2) $(NROBJ)
	$(BIND) -o $(EXE2) $(OBJ2) $(NROBJ) $(LIB)

$(EXE3) : $(OBJ3) $(IONOBJ)
	$(BIND) -o $(EXE3) $(OBJ3) $(IONOBJ) $(LIB)
#
# Define the dependencies for each object.
#
mglin.o : mglin.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) mglin.c

mgfas.o : mgfas.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) mgfas.c

ionosphere.o : ionosphere.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) ionosphere.c

IONOSPHERE_conductance.o : IONOSPHERE_conductance.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_conductance.c

IONOSPHERE_multigrid.o : IONOSPHERE_multigrid.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_multigrid.c

IONOSPHERE_fgrid.o : IONOSPHERE_fgrid.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_fgrid.c

IONOSPHERE_cgrid.o : IONOSPHERE_cgrid.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_cgrid.c

IONOSPHERE_restrict.o : IONOSPHERE_restrict.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_restrict.c

IONOSPHERE_prolong.o : IONOSPHERE_prolong.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_prolong.c

IONOSPHERE_smoother.o : IONOSPHERE_smoother.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_smoother.c

IONOSPHERE_residual.o : IONOSPHERE_residual.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_residual.c

IONOSPHERE_update.o : IONOSPHERE_update.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_update.c

IONOSPHERE_matcopy.o : IONOSPHERE_matcopy.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_matcopy.c

IONOSPHERE_matzero.o : IONOSPHERE_matzero.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_matzero.c

IONOSPHERE_util.o : IONOSPHERE_util.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_util.c

IONOSPHERE_magdim.o : IONOSPHERE_magdim.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_magdim.c

IONOSPHERE_magfac.o : IONOSPHERE_magfac.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_magfac.c

IONOSPHERE_magbfield.o : IONOSPHERE_magbfield.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_magbfield.c

IONOSPHERE_magvel.o : IONOSPHERE_magvel.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_magvel.c

IONOSPHERE_current.o : IONOSPHERE_current.c $(IONINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) IONOSPHERE_current.c

CPTG_util.o : CPTG_util.c $(CPTGINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) CPTG_util.c

NR_mglin.o : NR_mglin.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_mglin.c

NR_mgfas.o : NR_mgfas.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_mgfas.c

NR_addint.o : NR_addint.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_addint.c

NR_interp.o : NR_interp.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_interp.c

NR_lop.o : NR_lop.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_lop.c

NR_relax.o : NR_relax.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_relax.c

NR_relax2.o : NR_relax2.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_relax2.c

NR_resid.o : NR_resid.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_resid.c

NR_rstrct.o : NR_rstrct.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_rstrct.c

NR_slvsml.o : NR_slvsml.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_slvsml.c

NR_slvsm2.o : NR_slvsm2.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_slvsm2.c

NR_copy.o : NR_copy.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_copy.c

NR_matadd.o : NR_matadd.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_matadd.c

NR_matsub.o : NR_matsub.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_matsub.c

NR_fill0.o : NR_fill0.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_fill0.c

NR_anorm2.o : NR_anorm2.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_anorm2.c

NR_nrutil.o : NR_nrutil.c $(NRINC)
	$(CC) -c $(CPU) $(OPT) $(DEBUG) NR_nrutil.c
#
# Cleanup the multigrid solver home directory.
#
cleanup:
	rm -f *.o
	rm -f *~
