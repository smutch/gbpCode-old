############################################
#  Common Makefile for Greg Poole's code   #
#                                          #
#  Copyright (c) 2008  Gregory B. Poole,   #
#    Swinburne University of Technology,   #
#    Melbourne, Australia                  #
############################################

# Executables (use C99 standard so we have the trunc() function)
CC_NO_MPI  = gcc -std=c99
CC_NO_MPI  = icc 
CC_NO_MPI  = g++ -std=c++98
ifdef GBP_MPI
  CC_USE_MPI = $(GBP_MPI)/bin/mpic++
else
  CC_USE_MPI = mpic++
endif
MAKE       = make

# This ensures that we use standard (what is used in interactive shells) version of echo.
ECHO = /bin/echo

# This is needed to fix compilation errors re: undefined trunc() function
#CPPFLAGS := $(CPPFLAGS) -lm

# Compile getline() on Macs
ifndef USE_GETLINE
  USE_GETLINE=0
endif
ifneq ($(USE_GETLINE),0)
  CPPFLAGS := $(CPPFLAGS) -DUSE_GETLINE=1
else
  CPPFLAGS := $(CPPFLAGS) -DUSE_GETLINE=0
endif

# Set default so that real=float
ifndef USE_DOUBLE
  USE_DOUBLE=0
endif
ifneq ($(USE_DOUBLE),0)
  CPPFLAGS := $(CPPFLAGS) -DUSE_DOUBLE
endif
export USE_DOUBLE

# Set default MPI support
ifndef USE_MPI
  USE_MPI=0
endif
export USE_MPI
ifndef USE_MPI_IO
  USE_MPI_IO=0
endif
export USE_MPI_IO

# Set default Cuda support
ifndef USE_CUDA
  USE_CUDA=0
endif
export USE_CUDA

# Set default debugger support
ifndef USE_DEBUGGER
  USE_DEBUGGER=0
endif
export USE_DEBUGGER

# Set default file locations/destinations
ifndef GBP_SRC
  GBP_SRC:=$(PWD)
endif
export GBP_SRC
ifndef GBP_DAT
  GBP_DAT=$(GBP_SRC)/myData/
endif
export GBP_DAT
ifndef GBP_INC
  GBP_INC=$(GBP_SRC)/myInclude/
endif
export GBP_INC
ifndef GBP_LIB
  GBP_LIB=$(GBP_SRC)/myLib/
endif
export GBP_LIB
ifndef GBP_BIN
  GBP_BIN=$(GBP_SRC)/myBin/
endif
export GBP_BIN

# Are we using CUDA?
ifeq ($(USE_CUDA),1)
  CC_CUDA=nvcc
endif

# We need to tag-on the /cuda/ if USE_CUDA=1
ifeq ($(USE_CUDA),1)
  GBP_BIN_LOCAL:= $(GBP_BIN)/cuda/
  GBP_LIB_LOCAL:= $(GBP_LIB)/cuda/
else
  GBP_BIN_LOCAL:= $(GBP_BIN)/
  GBP_LIB_LOCAL:= $(GBP_LIB)/
endif

# We need to tag-on the /mpi/ if USE_MPI=1
ifeq ($(USE_MPI),1)
  GBP_BIN_LOCAL:= $(GBP_BIN_LOCAL)/mpi/
  GBP_LIB_LOCAL:= $(GBP_LIB_LOCAL)/mpi/
endif
export GBP_BIN_LOCAL
export GBP_LIB_LOCAL

# Set local variables
ifneq ($(wildcard Makefile.mylocal),)
  include Makefile.mylocal
else
  include Makefile.local
endif

# Set CPPFLAGS and LDFLAGS variables
CPPFLAGS := $(CPPFLAGS) -I$(GBP_INC) 
LDFLAGS := $(LDFLAGS) -L$(GBP_LIB_LOCAL)  

# Filter-out Cuda files if USE_CUDA!=1
ifneq ($(USE_CUDA),1)
	OBJEXCLUDE   :=$(patsubst %.o,%.ou,$(OBJFILES))
	OBJFILES_LIB1:=$(filter-out $(OBJEXCLUDE),$(OBJFILES))
else
	OBJFILES_LIB1:=$(OBJFILES)
endif
OBJFILES_LIB:=$(patsubst %.ou,%.o,$(OBJFILES_LIB1))

# Fill a file with a list of all the object files that
#   will contribute to this directory's archive (if defined)
ifneq ($(strip $(LIBFILE)),)
  LIBOBJSFILE=$(shell pwd)'/.'$(LIBFILE)'.objsfile'
  $(shell rm -f $(LIBOBJSFILE))
  $(shell touch $(LIBOBJSFILE))
  MAKE:=$(MAKE) LIBOBJSFILE=$(LIBOBJSFILE)
endif
LISTOBJ :=$(addprefix $(shell pwd)/,$(OBJFILES_LIB))' '
objfiles=$(shell cat $(LIBOBJSFILE))
#objfiles =$(call reverse,$(objfiles1))

# Set library information
export LIBS
include $(GBP_SRC)/Makefile.libs

# Turn on debugging information
ifeq ($(USE_DEBUGGER),1)
  CPPFLAGS := $(CPPFLAGS) -g
  LDFLAGS := $(LDFLAGS) -g
  ifeq ($(USE_CUDA),1)
    CPPFLAGS_CUDA := $(CPPFLAGS_CUDA) -G
  endif
else
  CPPFLAGS := $(CPPFLAGS) -O2
endif
export USE_DEBUGGER

#  Set compile flags, variables, etc. appropriately
LDFLAGS:=$(LDFLAGS) $(LIBS) 
ifeq ($(USE_MPI),1)
  CC           = $(CC_USE_MPI)
  INCTOUCH     = .install_inc
  DATATOUCH    = .install_data
  SCRIPTTOUCH  = .install_scripts
  LIBTOUCH1    = .install_lib_mpi
  BINTOUCH1    = .install_bin_mpi
  LIBTOUCH2    = .install_lib_nompi
  BINTOUCH2    = .install_bin_nompi
else
  CC           = $(CC_NO_MPI)
  INCTOUCH     = .install_inc
  DATATOUCH    = .install_data
  SCRIPTTOUCH  = .install_scripts
  LIBTOUCH1    = .install_lib_nompi
  BINTOUCH1    = .install_bin_nompi
  LIBTOUCH2    = .install_lib_mpi
  BINTOUCH2    = .install_bin_mpi
endif
CPPFLAGS := $(CPPFLAGS) -DGBP_DATA_DIR='"$(GBP_DAT)"'
OLDDATE=200701010101

################################
# TARGETS and RULES BEGIN HERE #
################################
.PHONY: clean force build_libobjsfile

# Common rules
all:     .print_status build_libobjsfile subdirs_all     $(INCTOUCH) $(LIBTOUCH1) $(BINTOUCH1) scripts_dir data_dir
include: .print_status build_libobjsfile subdirs_include $(INCTOUCH)
script:  .print_status build_libobjsfile subdirs_scripts scripts_dir
data:    .print_status build_libobjsfile subdirs_data    data_dir
lib:     .print_status build_libobjsfile subdirs_lib     $(LIBTOUCH1) 
bin:     .print_status build_libobjsfile subdirs_bin     $(BINTOUCH1)
clean:                 build_libobjsfile subdirs_clean
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
        	$(ECHO) -n "  " ; \
        	((i = i + 1)) ; \
	done
	@$(ECHO) -n "Cleaning-up..."
	@rm -rf .printed_status .print_status .install* .compile* .*.objsfile *.o *.ou *~ core.* *.a $(BINFILES) *.dSYM
	@$(ECHO) "Done."

build_libobjsfile:
ifdef LIBOBJSFILE
	@$(ECHO) -n $(LISTOBJ)" " >> $(LIBOBJSFILE)
endif

# Print a status message
.print_status: .printed_status
ifeq ($(MAKELEVEL),0)
	@$(ECHO)
	@$(ECHO) "Makefile flags:"
	@$(ECHO) "---------------"
ifneq ($(USE_DOUBLE),0)
	@$(ECHO) "USE_DOUBLE  is ON"
else
	@$(ECHO) "USE_DOUBLE  is OFF"
endif
ifneq ($(USE_MPI),0)
	@$(ECHO) "USE_MPI     is ON"
else
	@$(ECHO) "USE_MPI     is OFF"
endif
ifneq ($(USE_MPI_IO),0)
	@$(ECHO) "USE_MPI_IO  is ON"
else
	@$(ECHO) "USE_MPI_IO  is OFF"
endif
ifneq ($(USE_FFTW),0)
	@$(ECHO) "USE_FFTW    is ON"
else
	@$(ECHO) "USE_FFTW    is OFF"
endif
ifneq ($(USE_SPRNG),0)
	@$(ECHO) "USE_SPRNG   is ON"
else
	@$(ECHO) "USE_SPRNG   is OFF"
endif
ifneq ($(USE_HDF5),0)
	@$(ECHO) "USE_HDF5    is ON"
else
	@$(ECHO) "USE_HDF5    is OFF"
endif
ifneq ($(USE_CFITSIO),0)
	@$(ECHO) "USE_CFITSIO is ON"
else
	@$(ECHO) "USE_CFITSIO is OFF"
endif
ifneq ($(USE_CUDA),0)
	@$(ECHO) "USE_CUDA    is ON"
else
	@$(ECHO) "USE_CUDA    is OFF"
endif


        # Check if lib, include and bin directories exist.  If any do not, make them.
	@$(ECHO)
	@$(ECHO) "Directories:"
	@$(ECHO) "------------"
	@$(ECHO) -n "  Data      {"$(GBP_DAT)"} "
ifeq ($(wildcard $(GBP_DAT)),)
	@mkdir -p $(GBP_DAT)
	@$(ECHO) "(Created)"
else
	@$(ECHO) "(ok)"
endif
	@$(ECHO) -n "  Libraries {"$(GBP_LIB_LOCAL)"} "
ifeq ($(wildcard $(GBP_LIB_LOCAL)),)
	@mkdir -p $(GBP_LIB_LOCAL)
	@$(ECHO) "(Created)"
else
	@$(ECHO) "(ok)"
endif
	@$(ECHO) -n "  Headers   {"$(GBP_INC)"} "
ifeq ($(wildcard $(GBP_INC)),)
	@mkdir -p $(GBP_INC)
	@$(ECHO) "(Created)"
else
	@$(ECHO) "(ok)"
endif
	@$(ECHO) -n "  Binaries  {"$(GBP_BIN_LOCAL)"} "
ifeq ($(wildcard $(GBP_BIN_LOCAL)),)
	@mkdir -p $(GBP_BIN_LOCAL)
	@$(ECHO) "(Created)"
else
	@$(ECHO) "(ok)"
endif

	@$(ECHO)
	@$(ECHO) "Compile and linking parameters:"
	@$(ECHO) "-------------------------------"
	@$(ECHO) "CC     ="$(CC)
	@$(ECHO)
	@$(ECHO) "CPPFLAGS="$(CPPFLAGS)
	@$(ECHO)
	@$(ECHO) "LDFLAGS="'$(LDFLAGS)'
	@$(ECHO)
	@rm -rf .printed_status

endif
.printed_status:
	@touch  .printed_status

# Process subdirectories
.PHONY: subdirs_all subdirs_include subdirs_lib subdirs_bin subdirs_clean $(SUBDIRS)
subdirs_all:     $(addsuffix .GBPMAKEALL,    $(SUBDIRS))
subdirs_include: $(addsuffix .GBPMAKEINCLUDE,$(SUBDIRS))
subdirs_lib:     $(addsuffix .GBPMAKELIB,    $(SUBDIRS))
subdirs_bin:     $(addsuffix .GBPMAKEBIN,    $(SUBDIRS))
subdirs_clean:   $(addsuffix .GBPMAKECLEAN,  $(SUBDIRS))
$(addsuffix .GBPMAKEALL,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) all
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."
$(addsuffix .GBPMAKEINCLUDE,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) include
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."
$(addsuffix .GBPMAKESCRIPT,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) script
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."
$(addsuffix .GBPMAKEDATA,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) data
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."
$(addsuffix .GBPMAKELIB,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) lib
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."
$(addsuffix .GBPMAKEBIN,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) bin
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."
$(addsuffix .GBPMAKECLEAN,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) clean
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) "Done."

# Copy header files to the proper place
$(INCTOUCH): $(INCFILES)
ifneq ($(strip $(INCFILES)),)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) -n "Copying header file(s) {"$(INCFILES)"} to "$(GBP_INC)"..."
	@cp $(INCFILES) $(GBP_INC)
	@$(ECHO) "Done."
endif
	@touch $(INCTOUCH)

# Generate links to scripts
scripts_dir: $(addprefix $(GBP_BIN_LOCAL)/,$(SCRIPTS))
$(addprefix $(GBP_BIN_LOCAL)/,$(SCRIPTS)):
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) -n "Linking '"$(notdir $@)"' to bin directory..."
	@rm -rf $@
	@ln -s $(CURDIR)/$(notdir $@) $@
	@$(ECHO) "Done."

# Generate links to data files
data_dir: $(addprefix $(GBP_DAT)/,$(DATAFILES))
$(addprefix $(GBP_DAT)/,$(DATAFILES)):
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) -n "Linking '"$(notdir $@)"' to data directory..."
	@rm -rf $@
	@ln -s $(CURDIR)/data/$(notdir $@) $@
	@$(ECHO) "Done."

# Generate library
$(LIBTOUCH1): $(LIBFILE) $(OBJFILES_LIB1)
	@touch $(LIBTOUCH1)
$(LIBTOUCH2):
	@touch -t $(OLDDATE) $(LIBTOUCH2)
$(LIBFILE): $(OBJFILES_LIB1)
ifneq ($(strip $(LIBFILE)),)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ;  \
		((i = i + 1)) ; \
	done
	@$(ECHO) -n "Generating "$@"..."
	@ar ru $(LIBFILE) $(objfiles)
	@rm -f $(LIBOBJSFILE)
	@cp    $(LIBFILE) $(GBP_LIB_LOCAL)
endif

# Generate binaries
$(BINTOUCH1): $(BINFILES)
	@touch $(BINTOUCH1)
$(BINTOUCH2):
	@touch -t $(OLDDATE) $(BINTOUCH2)
$(BINFILES): $(LIBFILE) 
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) -n "Generating binary '"$@"'..."
ifneq ($(USE_CUDA),0)
	@$(CC_CUDA) $(CPPFLAGS) $(CPPFLAGS_CUDA) $@.c -o $@  $(LDFLAGS) -lcuda 
else
	@$(CC) $(CPPFLAGS) $@.c -o $@  $(LDFLAGS)
endif
	@cp $@ $(GBP_BIN_LOCAL)
	@$(ECHO) "Done."

# Generate object files (implicit rule)
%.o: %.c $(INCFILES) $(DEPENCIES) $(LIBTOUCH2)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
	@$(ECHO) -n "Compiling "$@"..."
	@$(CC) $(CPPFLAGS) -c $*.c
	@$(ECHO) "Done."

# Generate CUDA object files (implicit rule)
%.ou: %.cu $(INCFILES) $(DEPENCIES) $(LIBTOUCH2)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		$(ECHO) -n "  " ; \
		((i = i + 1)) ; \
	done
ifneq ($(USE_CUDA),0)
	@$(ECHO) -n "Compiling "$@"..."
	@$(CC_CUDA) $(CPPFLAGS) $(CPPFLAGS_CUDA) -c $*.cu
	@touch $*.ou
	@$(ECHO) "Done."
else
	 @$(ECHO) "Skipping "$@" (USE_CUDA is off)"
endif

