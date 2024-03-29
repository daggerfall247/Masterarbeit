################################################################################
#
# Makefile to compile and link C programs with or without MPI subroutines
#
# Version valid for Linux machines with MPICH
#
# "make" compiles and links the specified main programs and modules,
# using the specified libraries (if any), and produces the executables
#
# "make clean" removes all files generated by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all

# main programs and modules to be compiled

MAIN = meas_mf

RANDOM = ranlxd random_su3

ADMIN = init sun_vfunc endian error_checks utils

IO = IO_utils inp_IO config_IO

UPD = update metro exp_fct heatbath

OBS = plaquette wilson smearing 2pt topologicalcharge

DIRAC = dirac_wil spin_alg sources

INV = solv_cg

MODULES = $(RANDOM) $(ADMIN) $(IO) $(UPD) $(OBS) $(DIRAC) $(INV)


# search path for modules

MDIR = ../modules

VPATH = .:$(MDIR)/random:$(MDIR)/admin:$(MDIR)/io:$(MDIR)/update:$(MDIR)/obs:$(MDIR)/dirac:$(MDIR)/invert



# Logging option (-mpilog or -mpitrace or -mpianim)

LOGOPTION =


# additional include directories

INCPATH = ../include


# additional libraries

LIBS = m

LIBPATH = 


# scheduling and optimization options

# compiler flags
CFLAGS = -O3 -Wall

############################## do not change ###################################

SHELL=/bin/bash
CC=gcc
CLINKER=$(CC)

PGMS= $(MAIN) $(MODULES)

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(CC) -mcmodel=large $< -MM $(addprefix -I,$(INCPATH)) -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) -mcmodel=large $< -c $(CFLAGS) $(LOGOPTION) $(addprefix -I,$(INCPATH))


# rule to link object files

$(MAIN): %: %.o $(addsuffix .o,$(MODULES)) Makefile
	$(CLINKER) -mcmodel=large $< $(addsuffix .o,$(MODULES)) $(CFLAGS) $(LOGOPTION) \
        $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS)) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory

clean:
	@ -rm -rf *.d *.o *.alog *.clog *.slog $(MAIN)
.PHONY: clean

################################################################################
