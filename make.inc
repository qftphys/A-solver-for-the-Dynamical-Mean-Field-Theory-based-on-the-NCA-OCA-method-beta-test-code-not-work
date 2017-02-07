##$ COMPILER: suppprted compilers are ifort, gnu >v4.7
##$ use mpif90 for parallel compiler
FC=gfortran

##$ PLATFORM: supported platform are intel, gnu
##$ if using a parallel compiler (mpif90) check
##$ platform first using the cmd: $mpif90 -show
PLAT=gnu

##$ LOCATION OF THE scifor and dmft_tools DIRECTORIES
##$ is placed in different locations edit INCARGS here below
LIBDIR=/opt



##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
##$ if needed add yours to the list:
EXE=nca_test

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =.

##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin

ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz -prof-gen -prof-use
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debug extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp
endif

ifeq ($(PLAT),gnu)
INCARGS=-I$(LIBDIR)/scifor/gnu/include -L$(LIBDIR)/scifor/gnu/lib
INCARGS+=-I$(LIBDIR)/dmft_tools/gnu/include -L$(LIBDIR)/dmft_tools/gnu/lib
FFLAG = -O2 -funroll-all-loops -ffree-line-length-none $(INCARGS) 
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none  $(INCARGS)
#DFLAG = -O0 -p -g -Wall -fimplicit-none -ffpe-trap=zero,overflow,underflow  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -pedantic  -fbacktrace -fcheck=bounds -ffree-line-length-none  $(INCARGS)
OFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none $(INCARGS)
FPPFLAG =-cpp
endif



ARGS=  -ldmftt -lscifor
#-lfftpack -lminpack -llapack -lblas  -larpack -lparpack




##$ REVISION SOFTWARE VARIABLES
##$ 
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90
