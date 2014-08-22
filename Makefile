#COMPILER (PARALLEL)
FC=ifort

EXE=nca_test

DIR =./
DIREXE=./test

.SUFFIXES: .f90

# #REVISION SOFTWARE GIT:
# BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
# VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

OBJS= NCA_VARS_GLOBAL.o NCA_INPUT_VARS.o NCA_AUX_FUNX.o NCA_HLOCAL.o NCA_DIAG.o NCA_GREENS_FUNCTIONS.o NCA_MAIN.o DMFT_NCA.o

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

ARGS=-lscifor $(MKLARGS)


all:compile


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

.f90.o:	
	$(FC) $(FFLAG) -c $< 


completion:
	sf_lib_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

