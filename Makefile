EXEC   = darksage

SRCS   = 	./code/main.c \
		./code/core_read_parameter_file.c \
		./code/core_init.c \
		./code/core_io_tree.c \
		./code/core_cool_func.c \
		./code/core_build_model.c \
		./code/core_save.c \
		./code/core_mymalloc.c \
		./code/core_allvars.c \
		./code/model_infall.c \
		./code/model_cooling.c \
		./code/model_starformation_and_feedback.c \
		./code/model_disk_instability.c \
		./code/model_reincorporation.c \
		./code/model_mergers.c \
		./code/model_misc.c 

INCL   =	./code/core_allvars.h  \
			./code/core_proto.h  \
			./code/core_simulation.h  \
			./Makefile

OBJS  := $(SRCS:.c=.o)


# Neither of these options currently work!  Further development required.  Do not uncomment unless you intend on doing that development yourself.  See README for using MPI.
#USE-MPI=yes
#OPT += -DMINIMIZE_IO


UNAME := $(shell uname)

CC :=

## Set the C compiler if not set
ifeq ($(CC),)
  ## Make clang the default compiler on Mac
  ## But first check for clang-omp, use that if available
  ## clang/clang-omp is default on OSX
  ifeq ($(UNAME), Darwin)
    CC := clang
  else
    ## gcc is default on linux
    CC := gcc
  endif
endif

ifeq ($(CC),)
  $(error Error: Could not set compiler. Please either set "CC" in "Makefile" or via the command-line, "make CC=yourcompiler")
endif


ifdef USE-MPI
OPT += -DMPI  #  This creates an MPI version that can be used to process files in parallel
CC = mpicc  # sets the C-compiler
endif



OPTIMIZE = -O2 -Wall -Wextra -Wshadow   # optimization and warning flags (default)

# GSL automatic detection
GSL_FOUND := $(shell gsl-config --version 2>/dev/null)
ifndef GSL_FOUND
$(warning GSL not found in path - please install GSL before installing SAGE (or, update the PATH environment variable such that "gsl-config" is found))
# if the automatic detection fails, set GSL_DIR appropriately
GSL_DIR := /opt/local
GSL_INCL := -I$(GSL_DIR)/include
GSL_LIBDIR := $(GSL_DIR)/lib
# since GSL is not in PATH, the runtime environment might not be setup correctly either
# therefore, adding the compiletime library path is even more important (the -Xlinker bit)
GSL_LIBS := -L$(GSL_LIBDIR) -lgsl -lgslcblas -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
else
# GSL is probably configured correctly, pick up the locations automatically
GSL_INCL := $(shell gsl-config --cflags)
GSL_LIBDIR := $(shell gsl-config --prefix)/lib
GSL_LIBS   := $(shell gsl-config --libs) -Xlinker -rpath -Xlinker $(GSL_LIBDIR)
endif


LIBS   :=   -lm  $(GSL_LIBS)

CFLAGS :=   -g $(OPTIONS) $(OPT) $(OPTIMIZE) $(GSL_INCL)

default: all

$(EXEC): $(OBJS) 
#$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC) -g
	$(CC) $(OBJS) $(LIBS) -o $@

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) ./$(EXEC)

all:  tidy $(EXEC) clean
