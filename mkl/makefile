#===============================================================================
# Copyright 2004-2018 Intel Corporation.
#
# This software and the related documents are Intel copyrighted  materials,  and
# your use of  them is  governed by the  express license  under which  they were
# provided to you (License).  Unless the License provides otherwise, you may not
# use, modify, copy, publish, distribute,  disclose or transmit this software or
# the related documents without Intel's prior written permission.
#
# This software and the related documents  are provided as  is,  with no express
# or implied  warranties,  other  than those  that are  expressly stated  in the
# License.
#===============================================================================

##  Content:
##      Intel(R) Math Kernel Library Custom Shared Object builder
##
##******************************************************************************

help:
	@echo
	@echo "Custom Shared Object builder creates the universal MKL custom library"
	@echo
	@echo "Usage: make <target> [<options>]"
	@echo
	@echo "    target"
	@echo "      libuni   - static builder: uses static Intel(R) MKL interface+threading+core libraries"
	@echo "      help     - prints this help"
	@echo
	@echo "    <options>"
	@echo "        interface={lp64|ilp64} (for libuni only)"
	@echo "            Specifies programming interface for static builder only."
	@echo "            Default: lp64."
	@echo "        export=<file_name>"
	@echo "            The name of the file that contains the list of entry points"
	@echo "            to be included to the shared object."
	@echo "            Default: user_example_list (no extension)."
	@echo "        name=<dylib_name>"
	@echo "            The name of the shared object to be created. The extension .dylib will be added."
	@echo "            Default: mkl_custom."
	@echo "        xerbla=<err_handler>"
	@echo "            The name of the object file that contains the user's error handler."
	@echo "            By the default native Intel(R) MKL XERBLA is used."
	@echo "        threading={parallel|sequential}  (for static builder only)"
	@echo "            Specifies whether to use Intel(R) MKL in the threaded or sequential mode."
	@echo "            Default: parallel."
	@echo "        MKLROOT=<MKL_directory>"
	@echo "            Specifies the location of Intel(R) MKL libraries used to build the custom shared object."
	@echo "            Default: the Intel(R) MKL installation directory."
	@echo
	@echo "Usage examples:"
	@echo
	@echo "  make libuni"
	@echo "       Creates universal mkl_custom.dylib."
	@echo "       Function list is taken from the predefined file user_example_list."
	@echo "       The parallel Intel(R) MKL libraries are used to build mkl_custom.dylib."
	@echo "       Native Intel(R) MKL error handler is used."
	@echo

##------------------------------------------------------------------------------

ifndef MKLROOT
MKLROOT = ../..
endif

mkl_libpath="$(MKLROOT)/lib"
CMPLR_PATH="$(MKLROOT)/../compiler/lib"
BUILDER_DIR=$(MKLROOT)/tools/builder

#ifndef export
export=$(BUILDER_DIR)/user_example_list
#endif

#ifndef name
name=mkl_custom
#endif

ifndef interface
interface=lp64
endif

ifndef threading
threading=parallel
endif

IFACE_COMP_PART=intel
IFACE_THREADING_PART=intel

ifeq ($(interface),ilp64)
IFACE_LIB=libmkl_$(IFACE_COMP_PART)_ilp64.a
else
IFACE_LIB=libmkl_$(IFACE_COMP_PART)_lp64.a
endif

ifeq ($(threading),sequential)
THREADING_LIB=libmkl_sequential.a
OMP_LIB=
else
THREADING_LIB=libmkl_$(IFACE_THREADING_PART)_thread.a
OMP_LIB=-L$(CMPLR_PATH) -liomp5
endif
OMP_LIB += -lpthread

CORE_LIB=libmkl_core.a

PURE_IFACE_LIB=-L$(mkl_libpath) -lmkl_rt

LIBM=-lm

# ARCH_32=-arch i386
ARCH_64=-arch x86_64

fdir=$(BUILDER_DIR)/mkl_custom_fini
incs = -I$(MKLROOT)/include
mkl_custom_fini_32.o: $(fdir)/mkl_custom_fini.c ; gcc $(ARCH_32) $(incs) $^ -c -o $@
mkl_custom_fini_64.o: $(fdir)/mkl_custom_fini.c ; gcc $(ARCH_64) $(incs) $^ -c -o $@
mkl_custom_fini.o: $(if $(ARCH_32), mkl_custom_fini_32.o) $(if $(ARCH_64), mkl_custom_fini_64.o)
	lipo -create $^ -output $@

libuni ia32 intel64: check_export_file_name mkl_custom_fini.o
	@echo
	gcc $(ARCH_32) $(ARCH_64) -dynamiclib -single_module \
	$(xerbla) $(addprefix -u _,$(shell grep -v '[\#;]' $(export))) \
	mkl_custom_fini.o \
	$(mkl_libpath)/$(IFACE_LIB) \
	$(mkl_libpath)/$(THREADING_LIB) \
	$(mkl_libpath)/$(CORE_LIB) \
	$(OMP_LIB) $(LIBM) \
	-o "$(name).dylib" \
	-install_name "@rpath/$(name).dylib"
	@rm -f mkl_custom_fini.o mkl_custom_fini_32.o mkl_custom_fini_64.o

check_export_file_name: $(export)
$(export):
	@echo
	@echo ERROR: the export file '$(export) does not exist'
	@exit 1
