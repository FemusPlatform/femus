############################ MAKEFILE FOR FEMUS LIBRARY ##################################

# LIBRARY NAME ---------------------------------------------------------------------------
MAIN = $(FEMUS_DIR)/lib/libfemus_2d_$(METHOD).so


# COMPILER -------------------------------------------------------------------------------
CXX = mpicxx


# COMPILATION FLAGS ----------------------------------------------------------------------
ifeq ($(METHOD),opt)
  CFLAGS = -O3 -w -fPIC -std=c++11 -DNDEBUG
else
  CFLAGS = -O0 -felide-constructors -g -ansi -pedantic -W -Wall -Wextra \
  -Wno-long-long -Wpointer-arith -Wformat -Wparentheses -fPIC -std=c++11
endif


# INCLUDE DIRECTORIES --------------------------------------------------------------------
CONTRIB_HDF5_INCPATH =  $(HDF5_PATH)/include
FEMUS_CONTRIB        =  $(FEMUS_DIR)/contrib
MED_ROOT_DIR         =  $(MED_PATH)
med_ROOT_DIR         =  $(med_PATH)
MED_COUPL_DIR         =  $(MED_COUPL_PATH)

INCLUDES =  -I./
INCLUDES += -I$(FEMUS_DIR)/src/
INCLUDES += -I$(FEMUS_DIR)/solvers/ln_solvers/
INCLUDES += -I./DATA/
INCLUDES += -I$(FEMUS_DIR)/include/
INCLUDES += -I$(FEMUS_CONTRIB)/matrix/
INCLUDES += -I$(FEMUS_CONTRIB)/parallel/

INCLUDES += -I$(CONTRIB_HDF5_INCPATH)
INCLUDES += -I$(PETSC_DIR)/include
INCLUDES += -I${MED_ROOT_DIR}/include/salome/ 
INCLUDES += -I${med_ROOT_DIR}/include/
INCLUDES += -I${MED_COUPL_DIR}/include/


# SPECIAL INCLUDE FILES THAT REQUIRE RE-COMPILATION OF SOURCE FILES WHEN MODIFIED -------- 
INCLUDES_SRC  =  $(wildcard ./DATA/*.h)
INCLUDES_SRC +=  $(wildcard $(FEMUS_DIR)/include/*.h)
INCLUDES_CONTRIB  =  $(wildcard $(FEMUS_CONTRIB)/matrix/*.h)
INCLUDES_CONTRIB +=  $(wildcard $(FEMUS_CONTRIB)/parallel/*.h)


# LIBRARIES AND RELATIVE PATHS -----------------------------------------------------------
LFLAGS += -L$(HDF5_PATH)/lib/
LFLAGS += -L$(PETSC_DIR)/lib/
LFLAGS += -L$(MED_ROOT_DIR)/lib/salome/
LFLAGS += -L$(med_ROOT_DIR)/lib/
LFLAGS += -L$(MED_COUPL_DIR)/lib/

LIBS =  -lm  -lhdf5 -lpetsc -lmed -lmedloader -lmedcoupling \
	-linterpkernel -lmedcouplingremapper -lmedC 


# SOURCE FILES AND RELATIVE OBJECTS ------------------------------------------------------

# CONTRIB 
SRC_CONTRIB := $(sort $(notdir $(wildcard $(FEMUS_CONTRIB)/matrix/*.C)))
SRC_CONTRIB += $(sort $(notdir $(wildcard $(FEMUS_CONTRIB)/parallel/*.C)))
CONTRIB_OBJS = $(patsubst %.C, %-$(METHOD).o, $(SRC_CONTRIB))
LOC_CONTRIB_OBJS = $(addprefix obj/, $(CONTRIB_OBJS))

# FEMUS FILES
LIB_SRC_APP_FILES = $(sort $(notdir $(wildcard $(FEMUS_DIR)/src/*.C )))
LIB_SRC_APP_OBJS = $(patsubst %.C, %-$(METHOD).o, $(LIB_SRC_APP_FILES))
LIB_SRC_FEMUS_OBJS = $(addprefix obj/, $(LIB_SRC_APP_OBJS))


# MAKE COMMANDS --------------------------------------------------------------------------
# DISTINGUISH BETWEEN MAKE COMMANDS AND FILES
.PHONY: depend clean 

all:    $(MAIN)
	@echo  $(MAIN) succesfully compiled!

# RULE FOR LINKING OBJECT FILES
$(MAIN): $(LIB_SRC_FEMUS_OBJS) $(LOC_CONTRIB_OBJS)
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


# DEPENDENCY ON SPECIAL INCLUDES 
$(LIB_SRC_FEMUS_OBJS): $(INCLUDES_SRC)
$(LOC_CONTRIB_OBJS): $(INCLUDES_CONTRIB)


# RULE FOR LIBRARY FILES
obj/%-$(METHOD).o : %.C
	@mkdir -p $(dir $@)
	@echo "Compiling C++ (in "$(METHOD)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@		
	
vpath %.C $(FEMUS_DIR)/src/
vpath %.C $(FEMUS_CONTRIB)/matrix/
vpath %.C $(FEMUS_CONTRIB)/parallel/

# Remark 1 
# In this rule the automatic variables $< and $@ are used.
# The first is the name of the prerequisite (the .C file) while the $@ is the target of the
# rule (the .o file). By using these variables we can compile the two lists of sources 
# without defining anything else because we use an automatic suffix replacement. This rule 
# is very generic so be careful by using it.
#
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
DEPFILE  = .depend
DEPTOKEN = '\# MAKEDEPENDS'
DEPFLAGS = -Y -f $(DEPFILE) -s $(DEPTOKEN) -p ./


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
