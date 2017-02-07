include make.inc

OBJS= NCA_VARS_GLOBAL.o NCA_INPUT_VARS.o NCA_AUX_FUNX.o NCA_HLOCAL.o NCA_DIAG.o NCA_GREENS_FUNCTIONS.o NCA_MAIN.o DMFT_NCA.o

all: all version compile completion

all: FPPFLAG+=-D_


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FPPFLAG) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)$(BRANCH)

.f90.o:	
	$(FC) $(FPPFLAG) $(FFLAG) -c $< 

completion:
	scifor_completion.sh $(DIR)/$(EXE).f90

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

