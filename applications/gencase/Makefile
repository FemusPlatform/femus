################ MAKEFILE FOR FEMUS GENCASE ##################
# ----------------------------------------------------------
# --------------------- USES -------------------------------
# ----------------------------------------------------------
# make : compile the gencase
# make clean : delete the objects and the executable
#
#### -----------   FLAG -----------------------------------
# define the compiler to use
CXX = mpicxx
# define the compiling flags
ifeq ($(METHOD),opt)
  CFLAGS = -O3 -w -fPIC -std=c++14
else
  CFLAGS = -O0 -g -ansi -pedantic -W -Wall -Wextra -std=c++14 \
  -Wno-long-long -Wpointer-arith -Wformat -Wparentheses -fPIC 
endif


#### ------------   HDF5 -----------------------------
CONTRIB_HDF5_INCPATH=$(HDF5_PATH)/include
CONTRIB_HDF5_LIBPATH=$(HDF5_PATH)/lib
#### ------------   FEMUS -----------------------------
FEMUS_CONTRIB =$(FEMUS_DIR)/contrib
#### ------------   LIBMESH -----------------------------
LIBMESH_INCLUDE=$(LIBMESH_PATH)/include/
LIBMESH_CONTRIB=$(LIBMESH_PATH)/contrib/
LIBMESH_LIB=$(LIBMESH_PATH)/lib64/
#### ------------   PETSC -----------------------------
PETSC_PATH=$(PETSC_DIR)
PETSC_LIB=$(PETSC_PATH)/lib
#### ------------   MED -----------------------------
MED_ROOT_DIR         =  $(MED_PATH)
med_ROOT_DIR         =  $(med_PATH)
MED_COUPL_DIR        =  $(MED_COUPL_PATH)
#  --------------------  INCLUDE ----------------------------
# define any directories containing header files other than /usr/include
#
INCLUDES =  -I./
INCLUDES += -Isrc/
INCLUDES += -I$(APP_PATH)/DATA/
INCLUDES += -I$(FEMUS_DIR)/include/
INCLUDES += -I$(FEMUS_CONTRIB)/matrix/
INCLUDES += -I$(FEMUS_CONTRIB)/parallel/

INCLUDES += -I$(CONTRIB_HDF5_INCPATH)

INCLUDES += -I$(LIBMESH_INCLUDE)/base/  
INCLUDES += -I$(LIBMESH_INCLUDE)/ 
INCLUDES += -I$(LIBMESH_INCLUDE)/enums 
INCLUDES += -I$(LIBMESH_INCLUDE)/error_estimation
INCLUDES += -I$(LIBMESH_INCLUDE)/fe
INCLUDES += -I$(LIBMESH_INCLUDE)/libmesh
INCLUDES += -I$(LIBMESH_INCLUDE)/mesh
INCLUDES += -I$(LIBMESH_INCLUDE)/numerics
INCLUDES += -I$(LIBMESH_INCLUDE)/parallel
INCLUDES += -I$(LIBMESH_INCLUDE)/partitioning
INCLUDES += -I$(LIBMESH_INCLUDE)/physics
INCLUDES += -I$(LIBMESH_INCLUDE)/quadrature
INCLUDES += -I$(LIBMESH_INCLUDE)/reduced_basis
INCLUDES += -I$(LIBMESH_INCLUDE)/solution_transfer
INCLUDES += -I$(LIBMESH_INCLUDE)/solvers
INCLUDES += -I$(LIBMESH_INCLUDE)/geom
INCLUDES += -I$(LIBMESH_INCLUDE)/systems
INCLUDES += -I$(LIBMESH_INCLUDE)/utils

INCLUDES += -I$(PETSC_DIR)/include

INCLUDES += -I${MED_ROOT_DIR}/include/salome/ 
INCLUDES += -I${med_ROOT_DIR}/include/
INCLUDES += -I${MED_COUPL_DIR}/include/
# -------------------  LIBRARIES -------------------------------
# define library paths in addition to /usr/lib
LFLAGS := -L$(LIBMESH_LIB)
LFLAGS += -L$(CONTRIB_HDF5_LIBPATH)/
LFLAGS += -L$(FEMUS_DIR)
LFLAGS += -L$(PETSC_LIB)/
LFLAGS += -L$(LIBMESH_CONTRIB)/lib/x86_64-unknown-linux-gnu_opt/
LFLAGS += -L$(LIBMESH_CONTRIB)/tecplot/lib/x86_64-unknown-linux-gnu/

LFLAGS += -L$(MED_ROOT_DIR)/lib/salome/
LFLAGS += -L$(med_ROOT_DIR)/lib/
LFLAGS += -L$(MED_COUPL_DIR)/lib/
# and now specify the libraries' names
LIBS =  -lm  -lmesh_opt -lhdf5 -lpetsc -lmed -lmedloader -lmedcoupling \
	-linterpkernel -lmedcouplingremapper -lmedC 
# define any libraries to link into executable
# LIBS =  -lm -lmesh_opt -lhdf5 -lpetsc 
#-lHilbert  -lfparser  -lgzstream  -lmetis  -lnetcdf -lsfcurves \
#-ltriangle -lexodusii  -lgmv  -llaspack  -lnemesis  -lparmetis  -ltecio -ltetgen

# ----------------------------------------------------------
# ---------------- SOURCES and OBJECTS ---------------------
# ----------------------------------------------------------
# define the gencase source files and the include on which it depends
SRCS := $(wildcard ./src/*.C)
INCLUDES_SRCS :=  $(wildcard ./DATA/*.h)
# now define the correspondig objects of the source files simply by
# replacing the suffix .c of all words with the .o suffix.
# OBJS = $(patsubst %.C, %-$(METHOD).o, $(SRCS))

OBJ_DIR = ./obj/
SRC_FILES = $(sort $(notdir $(wildcard ./src/*.C )))
SRC_OBJS = $(patsubst %.C, %-$(METHOD).o, $(SRC_FILES))
OBJS_APP = $(addprefix $(OBJ_DIR), $(SRC_OBJS))


# finally define the executable file 
MAIN = $(FEMUS_DIR)/bin/gencase-$(METHOD)



# --------------------------------------------------------
# -------------- MAKE COMMANDS ---------------------------
# --------------------------------------------------------
# set PHONY variables to distinguish between file names 
# which could be in the directory and make rules
.PHONY: depend clean
# this is the executable file rule, the one executed when you call make
all:    $(SRC_FILES) $(MAIN)
	@echo   $(MAIN) succesfully compiled!


# this is the rule for linking the objects to obtain the executable file.
$(MAIN): $(OBJS_APP) 
	@echo $(SRCS)
	@echo $(SRC_FILES)
	@echo $(OBJS_APP)
	@echo "Linking "$@"..."
	@echo 
	@echo $^                | tr " " "\n"
	@echo $(LFLAGS)         | tr " " "\n"
	@echo $(LIBS)           | tr " " "\n"
	@echo "with compiling flags" $(CFLAGS)
	@echo "=========" $(MAIN) "==========="
	@echo 
	@$(CXX) $(CFLAGS) $(INCLUDES) -o $(MAIN) $^ $(LFLAGS) $(LIBS) 

# this is the rule for compiling the sources to obtain the objects.
# If the $(INCLUDES_SRCS) files change the gencase will be recompiled
%-$(METHOD).o : %.C %.h $(INCLUDES_SRCS)
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $<  -o  $@ 
	
# RULE FOR LOCAL FILES AND EQUATION MAP FILLER
obj/%.o : %.C
	@mkdir -p $(dir $@)
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@
	
obj/%-$(METHOD).o : %.C
	@mkdir -p $(dir $@)
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@
	
# vpath SETTINGS - WHERE TO LOOK FOR SOURCE FILES
vpath %.C ./src/


# Remark 1 
# In this rule the automatic variables $< and $@ are used.
# The first is the name of the prerequisite (the .C file) while the
# $@ is the target of the rule (the .o file). By using these variables
# we can compile the sources without defining anything else 
# because we use an automatic suffix replacement. This rule is very generic
# so be careful by using it.
# Remark 2 
# Writing a command in this way @command the command is not printed on the stdout.
# If you call @echo "Hello" it will be printed only 
# 
# Hello
# 
# if you omit the @ you get
# 
# echo "Hello"
# Hello
# 

clean:
	$(RM) ./src/*-$(METHOD).o *~ $(MAIN) *~ ./obj/*-$(METHOD).o
	


# -------------   Dependencies ------------------------
# include the dependency list 
DEPFILE	= .depend
DEPTOKEN	= '\# MAKEDEPENDS'
DEPFLAGS	= -Y -f $(DEPFILE) -s $(DEPTOKEN) -p ./


 depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)


 $(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CFLAGS) -- $(SRCS) >&/dev/null
	@echo "Updated .depend"




# "sinclude" allows for the Makefile to work when $(DEPFILE) is missing. 
sinclude  $(DEPFILE)

# DO NOT DELETE THIS LINE -- make depend needs it
