############################################
#  Common Makefile for Greg Poole's code   #
#                                          #
#  Copyright (c) 2008  Gregory B. Poole,   #
#    Swinburne University of Technology,   #
#    Melbourne, Australia                  #
############################################
# If you are running OSX-bash, you will    #
# want to fix a bug in echo.  You can do   #
# that with the following commands:        #
#   sudo port install bash                 #
#   sudo mv /bin/bash /bin/bash.old        #
#   sudo mv /bin/sh /bin/sh.old            #
#   sudo ln /opt/local/bin/bash /bin/bash  #
#   sudo ln /opt/local/bin/bash /bin/sh    #
############################################

# Executables
CC_NO_MPI  = gcc
ifndef GBP_MPI
  GBP_MPI=/usr/
endif
export GBP_MPI
CC_USE_MPI = $(GBP_MPI)/bin/mpicc
MAKE       = make

# Set default so that real=float
ifndef USE_DOUBLE
  USE_DOUBLE=0
endif
export USE_DOUBLE

# Set default to no MPI support
ifndef USE_MPI
  USE_MPI=0
endif
export USE_MPI

# Set default to no MPI-I/O support
ifndef USE_MPI_IO
  USE_MPI_IO=0
endif
export USE_MPI_IO

# Set default file locations/destinations
ifndef GBP_SRC
  GBP_SRC:=$(PWD)
endif
export GBP_SRC
ifndef GBP_DAT
  GBP_INC=$(GBP_DAT)/myData/
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

# We need to tag-on the /mpi/ if USE_MPI=1
ifeq ($(USE_MPI),1)
  GBP_BIN_LOCAL:= $(GBP_BIN)/mpi/
  GBP_LIB_LOCAL:= $(GBP_LIB)/mpi/
else
  GBP_BIN_LOCAL:= $(GBP_BIN)/
  GBP_LIB_LOCAL:= $(GBP_LIB)/
endif
export GBP_BIN_LOCAL
export GBP_LIB_LOCAL

# Set local variables
ifneq ($(wildcard Makefile.mylocal),)
  include Makefile.mylocal
else
  include Makefile.local
endif

# Start CCFLAGS and LDFLAGS variables
CCFLAGS = -I$(GBP_INC) -O2
LDFLAGS = -L$(GBP_LIB_LOCAL) 

# Fill a file with a list of all the object files that
#   will contribute to this directory's archive (if defined)
ifneq ($(strip $(LIBFILE)),)
  LIBOBJSFILE=$(shell pwd)'/.'$(LIBFILE)'.objsfile'
  $(shell rm -f $(LIBOBJSFILE))
  $(shell touch $(LIBOBJSFILE))
  MAKE:=$(MAKE) LIBOBJSFILE=$(LIBOBJSFILE)
endif
LISTOBJ:=$(addprefix $(shell pwd)/,$(OBJFILES))' '
objfiles = $(shell cat $(LIBOBJSFILE))

# Set library information
include $(GBP_SRC)/Makefile.libs

# Turn on debugging information
ifdef USE_DEBUGGER
  CCFLAGS := $(CCFLAGS) -g
  LDFLAGS := $(LDFLAGS) -g
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
  CCFLAGS     := $(CCFLAGS) -DUSE_MPI
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
ifneq ($(USE_MPI_IO),0)
  CCFLAGS  := $(CCFLAGS) -DUSE_MPI_IO
endif
ifneq ($(USE_DOUBLE),0)
  CCFLAGS := $(CCFLAGS) -DUSE_DOUBLE
endif
CCFLAGS := $(CCFLAGS) -DGBP_DATA_DIR='"$(GBP_DAT)"'
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
        	echo -n "  " ; \
        	((i = i + 1)) ; \
	done
	@echo -n "Cleaning-up..."
	@rm -rf .printed_status .print_status .install* .compile* .*.objsfile *.o *~ core.* *.a $(BINFILES)
	@echo "Done."

build_libobjsfile:
ifdef LIBOBJSFILE
	@echo -n $(LISTOBJ)" " >> $(LIBOBJSFILE)
endif

# Print a status message
.print_status: .printed_status
ifeq ($(MAKELEVEL),0)
	@echo
	@echo "Makefile flags:"
	@echo "---------------"
ifneq ($(USE_DOUBLE),0)
	@echo "USE_DOUBLE is ON"
else
	@echo "USE_DOUBLE is OFF"
endif
ifneq ($(USE_MPI),0)
	@echo "USE_MPI    is ON"
else
	@echo "USE_MPI    is OFF"
endif
ifneq ($(USE_MPI_IO),0)
	@echo "USE_MPI_IO is ON"
else
	@echo "USE_MPI_IO is OFF"
endif
ifneq ($(USE_FFTW),0)
	@echo "USE_FFTW   is ON"
else
	@echo "USE_FFTW   is OFF"
endif
ifneq ($(USE_SPRNG),0)
	@echo "USE_SPRNG  is ON"
else
	@echo "USE_SPRNG  is OFF"
endif

        # Check if lib, include and bin directories exist.  If any do not, make them.
	@echo
	@echo "Directories:"
	@echo "------------"
	@echo -n "  Data      {"$(GBP_DAT)"} "
ifeq ($(wildcard $(GBP_DAT)),)
	@mkdir -p $(GBP_DAT)
	@echo "(Created)"
else
	@echo "(ok)"
endif
	@echo -n "  Libraries {"$(GBP_LIB_LOCAL)"} "
ifeq ($(wildcard $(GBP_LIB_LOCAL)),)
	@mkdir -p $(GBP_LIB_LOCAL)
	@echo "(Created)"
else
	@echo "(ok)"
endif
	@echo -n "  Headers   {"$(GBP_INC)"} "
ifeq ($(wildcard $(GBP_INC)),)
	@mkdir -p $(GBP_INC)
	@echo "(Created)"
else
	@echo "(ok)"
endif
	@echo -n "  Binaries  {"$(GBP_BIN_LOCAL)"} "
ifeq ($(wildcard $(GBP_BIN_LOCAL)),)
	@mkdir -p $(GBP_BIN_LOCAL)
	@echo "(Created)"
else
	@echo "(ok)"
endif

	@echo
	@echo "Compile and linking parameters:"
	@echo "-------------------------------"
	@echo "CC     ="$(CC)
	@echo
	@echo "CCFLAGS="$(CCFLAGS)
	@echo
	@echo "LDFLAGS="'$(LDFLAGS)'
	@echo
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
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) all
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."
$(addsuffix .GBPMAKEINCLUDE,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) include
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."
$(addsuffix .GBPMAKESCRIPT,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) script
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."
$(addsuffix .GBPMAKEDATA,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) data
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."
$(addsuffix .GBPMAKELIB,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) lib
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."
$(addsuffix .GBPMAKEBIN,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) bin
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."
$(addsuffix .GBPMAKECLEAN,$(SUBDIRS)):
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Decending into subdirectory "$(basename $@)" ..."
	@$(MAKE) -s -C $(basename $@) clean
	@i=0 ; while [[ $$i -lt $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo "Done."

# Copy header files to the proper place
$(INCTOUCH): $(INCFILES)
ifneq ($(strip $(INCFILES)),)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo -n "Copying header file(s) {"$(INCFILES)"} to "$(GBP_INC)"..."
	@cp $(INCFILES) $(GBP_INC)
	@echo "Done."
endif
	@touch $(INCTOUCH)

# Generate links to scripts
scripts_dir: $(addprefix $(GBP_BIN_LOCAL)/,$(SCRIPTS))
$(addprefix $(GBP_BIN_LOCAL)/,$(SCRIPTS)):
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	echo -n "Linking '"$(notdir $@)"' to bin directory..."
	rm -rf $@
	ln -s $(CURDIR)/$(notdir $@) $@
	echo "Done."

# Generate links to data files
data_dir: $(addprefix $(GBP_DAT)/,$(DATAFILES))
$(addprefix $(GBP_DAT)/,$(DATAFILES)):
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo -n "Linking '"$(notdir $@)"' to data directory..."
	@rm -rf $@
	@ln -s $(CURDIR)/$(notdir $@) $@
	@echo "Done."

# Generate library
$(LIBTOUCH1): $(LIBFILE) $(OBJFILES)
	@touch $(LIBTOUCH1)
$(LIBTOUCH2):
	@touch -t $(OLDDATE) $(LIBTOUCH2)
$(LIBFILE): $(OBJFILES)
ifneq ($(strip $(LIBFILE)),)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		echo -n "  " ;  \
		((i = i + 1)) ; \
	done
	@echo -n "Generating "$@"..."
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
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo -n "Generating binary '"$@"'..."
	@$(CC) $(CCFLAGS) $@.c $(LDFLAGS) -o $@ 
	@cp $@ $(GBP_BIN_LOCAL)
	@echo "Done."

# Generate object files (implicit rule)
%.o: %.c $(INCFILES) $(DEPENCIES) $(LIBTOUCH2)
	@i=1 ; while [[ $$i -le $(MAKELEVEL) ]] ; do \
		echo -n "  " ; \
		((i = i + 1)) ; \
	done
	@echo -n "Compiling "$@"..."
	@$(CC) $(CCFLAGS) -c $*.c
	@echo "Done."

