#!/bin/csh
# IPE Makefile

# NOTE: Build parallelism can be controlled from the command line by setting the
#       MAKEJOBS variable, e.g.: MAKEJOBS=1 make theia_intel_parallel

# Static / Unconditional / Default variables
###nm20151113 DIRS          = main neutral eldyn flip plasma dummygptl driver
###nm20151110 DIRS          = main neutral flip eldyn plasma dummygptl driver
SHELL=/bin/csh


ifeq ($(TESTING),yes)
  DIRS        = dummygptl main integration_testing neutral eldyn flip plasma driver
else
  DIRS        = dummygptl main neutral eldyn flip plasma driver
endif

GPTL_INCFLAGS = -I../dummygptl/
GPTL_LDFLAGS  = -L../dummygptl -lgptl
PPP_FLAGS     = --Fmodule=module_decomp --Free --comment --HaloSize=10

TARGETS = clean                   \
          trillian_cray_parallel  \
          trillian_cray_serial    \
          trillian_intel_parallel \
          trillian_intel_serial   \
          trillian_pgi_parallel   \
          trillian_pgi_serial     \
          theia_intel_parallel    \
          theia_intel_serial \
          generic_pgi_serial \
          generic_gnu_serial \
          generic_intel_serial

# Dynamic / Conditional variables

ifeq ($(TESTING),yes)
  EXE         = $(CURDIR)/../bin/ipe_$(MACHINE)_$(COMPILER)_$(PARALLELISM)_testing.exe
else
  EXE         = $(CURDIR)/../bin/ipe_$(MACHINE)_$(COMPILER)_$(PARALLELISM).exe
endif
FFLAGS        = $(strip $(OPT_FLAGS) $(GPTL_INCFLAGS) $(GPTL_LDFLAGS))
PARALLEL_MAKE = $(MAKE) -j $(MAKEJOBS)
PPP           = $(SMS)/bin/ppp
SMS_INCFLAGS  = -I$(SMS)/include
###nm20151014 SMS_LDFLAGS   = -L$(SMS)/lib -lsms
SMS_LDFLAGS   = -L$(SMS)/lib -lsms  $(GPTL_LDFLAGS)

ifndef DEBUG
  DEBUG = no
endif

# generic stanza

ifeq ($(MACHINE),generic)
  GPTL_USE_GPTL = yes
ifndef MAKEJOBS
  MAKEJOBS      = 1
endif

# PGI
  ifeq ($(COMPILER),pgi)
    MODCMD    = $(shell ./modcmd generic pgi)
    FCP       = pgfortran
    FCS       = pgfortran
    OPT_FLAGS = -fast -Mpreprocess
    SMS       =
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -Mbounds -Mpreprocess
    endif

# GNU
  else ifeq ($(COMPILER),gnu)
    MODCMD    = $(shell ./modcmd generic pgi)
    FCP       = gfortran -ffree-line-length-none
    FCS       = gfortran -ffree-line-length-none
    OPT_FLAGS = -Ofast -cpp
    SMS       =
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -Mbounds -cpp
    endif

#INTEL
  else ifeq ($(COMPILER),intel)
    MODCMD    = $(shell ./modcmd generic pgi)
    FCP       = ifort
    FCS       = ifort
    #OPT_FLAGS = -O2 -qopt-report-phase=vec -qopt-report=2
    OPT_FLAGS = -fast -fpp
    SMS       =
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -Mbounds -fpp
    endif
  endif
endif

# trillian stanza

ifeq ($(MACHINE),trillian)
  GPTL_USE_GPTL = no
ifndef MAKEJOBS
  MAKEJOBS      = 12
endif
  ifeq ($(GPTL_USE_GPTL),yes)
    $(error GPTL support currently unavailable on trillian)
  endif
  ifeq ($(COMPILER),cray)
    MODCMD    = $(shell ./modcmd trillian cray)
    FCP       = ftn
    FCS       = ftn
    OPT_FLAGS = -O3 -em -ef -eZ -ra -h fp0
    SMS       = /mnt/lustre/lus0/space/madden/SMS/legacy_r689/trillian-cray
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -em -ef -eZ -ra
    endif
  endif
  ifeq ($(COMPILER),intel)
    MODCMD    = $(shell ./modcmd trillian intel)
    FCP       = ftn
    FCS       = ftn
    OPT_FLAGS = -O3 -fp-model precise -fpp
    SMS       = /mnt/lustre/lus0/space/madden/SMS/legacy_r689/trillian-intel
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -ftz -traceback -check bounds -fpe0 -fpp
    endif
  endif
  ifeq ($(COMPILER),pgi)
    MODCMD    = $(shell ./modcmd trillian pgi)
    FCP       = ftn
    FCS       = ftn
    OPT_FLAGS = -O2 -Mpreprocess
    SMS       = /mnt/lustre/lus0/space/madden/SMS/legacy_r689/trillian-pgi
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -Mbounds -Mpreprocess
    endif
  endif
endif

# theia stanza

ifeq ($(MACHINE),theia)
  GPTL_USE_GPTL = yes
ifndef MAKEJOBS
  MAKEJOBS      = 24
endif 
  ifeq ($(COMPILER),intel)
    MODCMD    = $(shell ./modcmd theia intel)
    FCP       = mpiifort
    FCS       = ifort
    CC        = icc
    OPT_FLAGS = -O3 -fp-model precise  -vec-report2 -fpp 
    SMS       = /contrib/sms/legacy/SMSr856-intel
    ifeq ($(DEBUG),yes)
      OPT_FLAGS = -g -O0 -ftz -traceback -check bounds -fpe0 -fpp
    endif
    ifeq ($(PARALLELISM),parallel)
      ifeq ($(GPTL_USE_GPTL),yes)
        GPTL_INCFLAGS = -I/contrib/gptl/gptl-v5.4.4_impi_noomp/include
        GPTL_LDFLAGS  = -L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl
      endif
    endif
    ifeq ($(PARALLELISM),serial)
      ifeq ($(GPTL_USE_GPTL),yes)
#       $(error GPTL support currently unavailable on theia serial intel)
        GPTL_INCFLAGS = -I/contrib/gptl/gptl-v5.4.4_impi_noomp/include
        GPTL_LDFLAGS  = -L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl
      endif
    endif
  endif
endif


# This statement is machine/compiler independent
ifeq ($(TESTING),yes)
  OPT_FLAGS+= -DTESTING
endif

VARS = SHELL="$(SHELL)"               \
       DIRS="$(DIRS)"                 \
       EXE="$(EXE)"                   \
       FCP="$(FCP)"                   \
       FCS="$(FCS)"                   \
       FFLAGS="$(FFLAGS)"             \
       MODCMD="$(MODCMD)"             \
       PARALLELISM="$(PARALLELISM)"   \
       PPP="$(PPP)"                   \
       PPP_FLAGS="$(PPP_FLAGS)"       \
       SMS_INCFLAGS="$(SMS_INCFLAGS)" \
       SMS_LDFLAGS="$(SMS_LDFLAGS)"   \
       GPTL_LDFLAGS="$(GPTL_LDFLAGS)"

# Targets

.PHONY: $(TARGETS) all ipe

all:
	$(info Valid targets are:)
	$(foreach x,$(TARGETS),$(info - $(x)))
	$(error Invalid target specified)

trillian_cray_parallel:
	$(MAKE) ipe MACHINE=trillian COMPILER=cray PARALLELISM=parallel

trillian_cray_serial:
	$(MAKE) ipe MACHINE=trillian COMPILER=cray PARALLELISM=serial

trillian_intel_parallel:
	$(MAKE) ipe MACHINE=trillian COMPILER=intel PARALLELISM=parallel

trillian_intel_serial:
	$(MAKE) ipe MACHINE=trillian COMPILER=intel PARALLELISM=serial

trillian_pgi_parallel:
	$(MAKE) ipe MACHINE=trillian COMPILER=pgi PARALLELISM=parallel

trillian_pgi_serial:
	$(MAKE) ipe MACHINE=trillian COMPILER=pgi PARALLELISM=serial

theia_intel_parallel:
ifeq ($(TESTING),yes)
	$(MAKE) cmdi MACHINE=generic COMPILER=intel PARALLELISM=serial
	make --directory=integration_testing/ clean
endif
	$(MAKE) ipe MACHINE=theia COMPILER=intel PARALLELISM=parallel

theia_intel_serial:
ifeq ($(TESTING),yes)
	$(MAKE) cmdi MACHINE=generic COMPILER=intel PARALLELISM=serial
	make --directory=integration_testing/ clean
endif
	$(MAKE) ipe MACHINE=theia COMPILER=intel PARALLELISM=serial

generic_pgi_serial:
ifeq ($(TESTING),yes)
	$(MAKE) cmdi MACHINE=generic COMPILER=pgi PARALLELISM=serial
	make --directory=integration_testing/ clean
endif
	$(MAKE) ipe MACHINE=generic COMPILER=pgi PARALLELISM=serial

generic_gnu_serial:
ifeq ($(TESTING),yes)
	$(MAKE) cmdi MACHINE=generic COMPILER=gnu PARALLELISM=serial
	make --directory=integration_testing/ clean
endif
	$(MAKE) ipe MACHINE=generic COMPILER=gnu PARALLELISM=serial

generic_intel_serial:
ifeq ($(TESTING),yes)
	$(MAKE) cmdi MACHINE=generic COMPILER=intel PARALLELISM=serial
	make --directory=integration_testing/ clean
endif
	$(MAKE) ipe MACHINE=generic COMPILER=intel PARALLELISM=serial

ipe:
	$(if $(filter 0,$(MAKELEVEL)),$(call USAGE))
	$(foreach x,$(DIRS),$(PARALLEL_MAKE) -C $(x) $(VARS) &&) true

cmdi:
	make --directory=integration_testing/ -f cmdi.mak cmdi $(VARS)
	mv integration_testing/cmdi ../bin/

clean:
	$(foreach x,$(DIRS),$(MAKE) -C $(x) clean &&) true
	make --directory=integration_testing/ clean
