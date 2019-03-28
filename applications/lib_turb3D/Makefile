######### MAKEFILE FOR FEMUS APPLICATION ##################
# ----------------------------------------------------------
# --------------------- USES -------------------------------
# ----------------------------------------------------------
# make : compile the main application we are running
# make gencase : compile the gencase
# make clean : delete the objects files and the executable just in the app directory
# make src_clean : delete everything compiled in femus
# make resu_clean : delete the solutions in the app RESU dir, excluding mesh and matrixes
# make all_clean : delete everything compiled in femus and clean the app RESU directory


# ----------------------------------------------------------
# -------------------- FLAGS -------------------------------
# ----------------------------------------------------------

OBJ_DIR = obj/
TURB_SOLVERS = $(FEMUS_DIR)/solvers/MGSolverTURB/solvers/
TURB_UTILITY = $(FEMUS_DIR)/solvers/MGSolverTURB/utility/
TURB_INCLUDE = $(FEMUS_DIR)/solvers/MGSolverTURB/include/


CXX = mpicxx


ifeq ($(METHOD),opt)
  CFLAGS = -O3 -w -fPIC -std=c++11 -DNDEBUG
else
  CFLAGS = -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra \
  -Wno-long-long -Wpointer-arith -Wformat -Wparentheses -fPIC -std=c++11
endif


# ----------------------------------------------------------------------------------------
# PATHS FOR HEADER FILES
# ----------------------------------------------------------------------------------------
CONTRIB_HDF5_INCPATH =  $(HDF5_PATH)/include
FEMUS_CONTRIB        =  $(FEMUS_DIR)/contrib
MED_ROOT_DIR         =  $(MED_PATH)
med_ROOT_DIR         =  $(med_PATH)
MED_COUPL_DIR         =  $(MED_COUPL_PATH)

INCLUDES =  -I$(TURB_INCLUDE)
INCLUDES += -I./DATA/
INCLUDES += -I$(FEMUS_DIR)/include/
INCLUDES += -I$(FEMUS_CONTRIB)/matrix/
INCLUDES += -I$(FEMUS_CONTRIB)/parallel/
INCLUDES += -I$(CONTRIB_HDF5_INCPATH)
INCLUDES += -I$(PETSC_DIR)/include
INCLUDES += -I${MED_ROOT_DIR}/include/salome/ 
INCLUDES += -I${med_ROOT_DIR}/include/
INCLUDES += -I${MED_COUPL_DIR}/include/


# ----------------------------------------------------------------------------------------
# LIST OF HEADER FILES THAT REQUIRE LIBRARY COMPILATION WHEN MODIFIED
# ----------------------------------------------------------------------------------------
INCLUDES_SRC  =  $(wildcard ./DATA/*.h)
INCLUDES_SRC +=  $(wildcard $(FEMUS_DIR)/include/*.h)
INCLUDES_SRC +=  $(wildcard $(TURB_INCLUDE)/*.h)
INCLUDES_CONTRIB  =  $(wildcard $(FEMUS_CONTRIB)/matrix/*.h)
INCLUDES_CONTRIB +=  $(wildcard $(FEMUS_CONTRIB)/parallel/*.h)


# ----------------------------------------------------------------------------------------
# PATHS TO NEEDED LIBRARIES
# ----------------------------------------------------------------------------------------
LFLAGS += -L$(HDF5_PATH)/lib/
LFLAGS += -L$(PETSC_DIR)/lib/
LFLAGS += -L$(MED_ROOT_DIR)/lib/salome/
LFLAGS += -L$(med_ROOT_DIR)/lib/
LFLAGS += -L$(MED_COUPL_DIR)/lib/


# ----------------------------------------------------------------------------------------
# NEEDED LIBRARIES
# ----------------------------------------------------------------------------------------
LIBS =  -lm  -lhdf5 -lpetsc -lmed -lmedloader -lmedcoupling \
	-linterpkernel -lmedcouplingremapper -lmedC 


# ----------------------------------------------------------------------------------------
# SOURCE FILES TO BE COMPILED
# ----------------------------------------------------------------------------------------
SRC_APP_FILES = $(sort $(notdir $(wildcard $(TURB_SOLVERS)/*.C )))
SRC_APP_FILES += $(sort $(notdir $(wildcard $(TURB_UTILITY)/*.C )))
SRC_APP_OBJS = $(patsubst %.C, %-$(METHOD).o, $(SRC_APP_FILES))
OBJS_APP = $(addprefix $(OBJ_DIR), $(SRC_APP_OBJS))

MAIN = $(FEMUS_DIR)/lib/libturb_3d_$(METHOD).so


# --------------------------------------------------------
# -------------- MAKE COMMANDS ---------------------------
# --------------------------------------------------------
# set PHONY variables to distinguish between file names 
# which could be in the directory and make rules
.PHONY: depend clean 
# this is the executable file rule, the one executed when you call make
all:    $(MAIN)
	@echo  $(MAIN) succesfully compiled!

# this is the rule for linking the objects to obtain the executable file.
# There are two groups of objects, the SRC and the CONTRIB objects, which
# are recompiled in different cases
$(MAIN): $(OBJS_APP) 
	@echo 
	@echo "Linking "$@"..."
	@echo 
	@echo $^          | tr " " "\n"
	@echo $(LFLAGS)   | tr " " "\n"
	@echo $(LIBS)     | tr " " "\n"
	@echo "with compiling flags" $(CFLAGS)
	@echo "=========" $(MAIN) "==========="
	@echo 
	@$(CXX) -shared $(CFLAGS) $(INCLUDES) -o $(MAIN) $^ $(LFLAGS) $(LIBS)

# this is the rule for compiling the sources to obtain the objects.
# The actual rule is the third line below, which is very generic,
# while the first two lines define the dependencies of the two different
# groups of objects. Note that the SRC_OBJS are recompiled if the
# INCLUDES_SRC files change, while the CONTRIB_OBJS are recompiled if the
# INCLUDES_CONTRIB files change.
$(OBJS_APP): $(INCLUDES_SRC)
%-$(METHOD).o : %.C
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

# RULE FOR LIBRARY FILES
obj/%-$(METHOD).o : %.C
	@mkdir -p $(dir $@)
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@
	


vpath %.C $(FEMUS_DIR)/solvers/MGSolverTURB/solvers/
vpath %.C $(FEMUS_DIR)/solvers/MGSolverTURB/utility/




# Remark 1 
# In this rule the automatic variables $< and $@ are used.
# The first is the name of the prerequisite (the .C file) while the
# $@ is the target of the rule (the .o file). By using these variables
# we can compile the two lists of sources without defining anything else 
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

# The other useful rules below are explained in the MAKEFILE uses, see above
clean:
	$(RM) obj/*-$(METHOD).o *~ $(MAIN)

# -------------   MAKE depend and depend file ------------------------
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
