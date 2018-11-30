# ===========================================================================================
# =================  CONFIGURATION SCRIPT FOR FEMus APPLICATIONs  1/10/2017 =================
# ===========================================================================================
#!/bin/sh
# In order to use femus run source configure.sh from femus directory.
# 
# $$ source configure.sh APPNAME_number options
# 
# APPNAME must be an existing application in the template_appl directory,
# else you can create an empty application to write on your own.
# Open always three tabs for your application, one in dbg mode, one 
# in opt mode and the last with the correct gencase.
# ===========================================================================================
########################################################################
########## Over here there MUST be the definition of       #############
########## $INSTALLATION_DIR, given when installing femus; #############
########## otherwise, write it on your own                 #############
########################################################################

# function to run FEMuS 
# the function takes as input the number of processors
# a log file is printed by redirecting the output of std::cerr
function runFEMuS {
   make -j$1 
   mpiexec -np $1 $FM_MYAPP-opt 2> messages.log
}


# ===========================================================================================
#     Environment platform settings: $PATH  +  $LD_LIBRARY_PATH
# ===========================================================================================
function CreateEnvironment {
 
  export MYHOME=$PWD
  
  # set $1= INSTALLATION_DIR=$PWD
  echo "INSTALLAION_PATH=" $1
   # ----------------------------------------------------------------------------------------
   #                         PATH for LIBRARIES  
   # ----------------------------------------------------------------------------------------
   export FEMUS_DIR=$1$FEMUS_NAME"/"
#    export FEMUS_DIR=$PWD
   echo "FEMUS_DIR - create environment: " $FEMUS_DIR
   cd $FEMUS_DIR
   # -----------  LIBMESH ---------------
   export LIBMESH_PATH=$1/libmesh/
   echo "LIBMESH_PATH" $LIBMESH_PATH
   #  -----------  PETSC   ---------------
   export PETSC_DIR=$1/petsc/
  echo "PETSC_DIR" $PETSC_DIR
  # -----------  OPENMPI  ---------------
  export MPI_BIN_PATH=$1/openmpi/bin/
  export MPI_LIB_PATH=$1/openmpi/lib64/
  echo "MPI_LIB_PAT" $MPI_LIB_PATH
  echo "MPI_BIN_PATH" $MPI_BIN_PATH
  # -----------  HDF5  ---------------
  export HDF5_PATH=$1/hdf5/
  echo "HDF5_PATH" $HDF5_PATH
  # -----------  med  ---------------
  export salome_dir=$1/salome/
  export PATH=$salome_dir:$PATH
  export med_PATH=$1/salome/med
  echo "med_PATH" $med_PATH
  # -----------  MED  ---------------
  export MED_PATH=$1/salome/MED_mod
  echo "MED_PATH" $MED_PATH
  export MED_COUPL_PATH=$INSTALLATION_DIR/salome/MED_coupling
  echo "MED_COUPL_PATH" $MED_COUPL_PATH
  # -----------------------------------------------------------------------------------------
  #                             GUI 
  # -----------------------------------------------------------------------------------------
  export FM_GUI=gui
  GUI_FEMUS="../../"$FM_GUI
  GUI_CMD="code_femus"
  # -----------------------------------------------------------------------------------------
  #                         LIBRARIES 
  # -----------------------------------------------------------------------------------------
  # ------- MPI library ---------------------------------------------------------------------
  export PATH=$MPI_BIN_PATH:$PATH
  export LD_LIBRARY_PATH=$MPI_LIB_PATH:$LD_LIBRARY_PATH
  # ------- Libmesh library -----------------------------------------------------------------
  export LD_LIBRARY_PATH=$LIBMESH_PATH/lib64/:$LD_LIBRARY_PATH
  # ------- PETSC library -------------------------------------------------------------------
  export LD_LIBRARY_PATH=$PETSC_DIR/lib/:$LD_LIBRARY_PATH
  # ------- HDF5 library --------------------------------------------------------------------
  export LD_LIBRARY_PATH=$HDF5_PATH/lib:$LD_LIBRARY_PATH
  # ------- med tool library ----------------------------------------------------------------
  export LD_LIBRARY_PATH=$med_PATH/lib:$LD_LIBRARY_PATH
  # ------- MED module library --------------------------------------------------------------
  export LD_LIBRARY_PATH=$MED_PATH/lib/salome:$LD_LIBRARY_PATH
  # ------- Med coupling library ------------------------------------------------------------
  export LD_LIBRARY_PATH=$MED_COUPL_PATH/lib:$LD_LIBRARY_PATH
  # ------- FEMUS ---------------------------------------------------------------------------
  export LD_LIBRARY_PATH=$FEMUS_DIR/contrib/laspack:$LD_LIBRARY_PATH  #only for VOF
  export LD_LIBRARY_PATH=$FEMUS_DIR/lib:$LD_LIBRARY_PATH  #only for VOF
}

 
  
# ===========================================================================================================
# function CheckExistence: check first argument $1 as name application
# ===========================================================================================================
function CheckExistence {
#     exist or  create  USER_APPL/$FM_MYAPP directory
   echo "current= " $PWD
   cd $FEMUS_DIR
   if [ -d  USER_APPL/$1 ]; then
     # USER_APPL/$FM_MYAPP  exists
     echo "USER_APPL/"$1 " directory exists  " 
   else
     # USER_APPL/$FM_MYAPP does not exist
     echo "USER_APPL/"$1 " directory does not  exist  "
     echo "Do you want to create the USER_APPL/"$1  "application ? [yes/no]"
     read test_app
     if test "$test_app" = "yes"; then 
       cp -r template_appl/$APP_NAME/   USER_APPL/$1
     else
       return
     fi
   fi
     cd USER_APPL/$1
     echo "APPLICATION=" $1 "APP_NAME=" $APP_NAME "APP_VER=" $APP_VER "METHOD=" $METHOD "in USER_APPL/"$FM_MYAPP
     echo  "path " $PWD 
     export APP_PATH=$PWD 
}

# ===========================================================================================================
# function CreateApp: create application with name $1 with mode dbg,opt, etc....
# ===========================================================================================================
function CreateApp {
if test "$1" = "makelib"; then
   cd $FEMUS_DIR/lib/lib_femus
   return;
else
  CheckExistence $1
  if test "$2" = "dbg"; then
   export METHOD=dbg
   OPTION="nogui"
   export PETSC_ARCH=linux-dbg
  else
   export METHOD=opt
   export PETSC_ARCH=linux-opt   
  fi      
  ## Rename the tab as APP_NAME + METHOD
  echo $'\033]30;' $1 $2 $'\007' 
  
  if test "$OPTION" = "nogui"; then
   return;
  fi
  if test "$OPTION" = "gencase"; then
   $PWD
   cd ../gencase
   return;
  fi
  if test "$OPTION" = "interpolator"; then
   $PWD
   cd ../Interpolator
   printf "\033c"
   echo "
   
   ${red}${bold} In order to use the Interpolation application you need to define 
          \"HAVE_MED\" 
 inside the file" $FM_MYAPP"/DATA/Solverlib_conf.h 
 -) Make a resu_clean
 -) Set restart to 0
 -) Set the new meshes in \"param_files\"
 -) Run gencase with n processors
 -) Run Interpolator with n processors
 -) Set restart to 1000 
 -) Restart your application${reset}
 
 
 "
   return;
  fi
fi  
}

# ===========================================================================================================
# function ApplExistence: check first argument $1 as application name in template directory
# ===========================================================================================================
function ApplExistence {
APP=$1
APP_NAME=${APP/%_*/""}
APP_VER=${APP/#*_/""}
echo  "APP_NAME" $APP_NAME "APP_VER" $APP_VER
if test "$APP_NAME" = "$1";  then APP_NAME="None";fi
echo  "APP_NAME" $APP_NAME "APP_VER" $APP_VER
cd template_appl 
ls
if [ -d  $APP_NAME ]; then
  echo $APP_NAME " is an application template - ok  "
  cd ..
else
  echo  " "$1" name is invalid: not VALID application template or _n is missing "
  cd ..
    echo
    echo " Use: source configure.sh application_n  option"
    echo
    echo
    echo
    echo "   applications: applications_n" 
    echo "                 example nse_1, fsi_2, ..; n=version "
    echo
    echo
    echo "   option      : None, opt, dbg, gui, nogui, gencase "
    echo "		   None         = no gui in opt mode "
    echo "		   opt          = salome gui in opt mode"
    echo "		   dbg          = no gui in dbg mode"
    echo "		   nogui        = no gui in opt mode"
    echo "		   gencase      = no gui opt mode in gencase directory"
    echo "                   interpolator = application for interpolation of solution fields"
    echo "                   all          = open 4 tabs: dbg, nogui, gencase and interpolator"
    echo
fi
  }
  
  
# ===========================================================================================================
# main script
# ===========================================================================================================
# Argument script: $1 = appilcation name    $2=option 
# ===========================================================================================================
  red=`tput setaf 1`
  bold=`tput bold `
  green=`tput setaf 2`
  reset=`tput sgr0`
# ===========================================================================================================
# ===========================================================================================================
  # installation path=software_dir
  export FEMUS_PATH=$PWD
  
  cd ../

  export INSTALLATION_DIR=$PWD
  
  export LL=`expr length $INSTALLATION_DIR`
  export LLL=`expr length $FEMUS_PATH`
  export L1=`expr $LLL - $LL`
  
  export FEMUS_NAME=`expr substr $FEMUS_PATH $(($LL +1)) $L1`

  echo "---------------------------------------------------------------"
  echo "FEMUS PATH   " $FEMUS_PATH
  echo "INSTALLATION DIR   "$INSTALLATION_DIR
  echo "FEMUS NAME   " $FEMUS_NAME
  echo "---------------------------------------------------------------"
  
  echo "INSTALLAION_PATH=" $INSTALLATION_DIR
  # environment MPI PETSC libmesh HDF5 MED Medcoupling ......
  CreateEnvironment  $INSTALLATION_DIR
  export curdir=$PWD
  # Application with name $1: check (CheckExistence) and creation (CreateApp)
  # ApplExistence -> CheckExistence -> CreateApp ......
  ApplExistence $1
  export FM_MYAPP=$1
  export OPTION=$2 

 
 export CASE="nocase"
 
  if test "$OPTION" = "all"; then
    # multiple tabs
    var1=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession) 
#     echo $var1
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var1)) runCommand "cd  $curdir && source configure.sh $FM_MYAPP nogui"
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
    var2=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession) 
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var2)) runCommand "cd $curdir  && source configure.sh $FM_MYAPP gencase"
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
    var3=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession)   
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var3)) runCommand "cd $curdir  && source configure.sh $FM_MYAPP interpolator"
    var4=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession) 
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) runCommand "cd $curdir  && source configure.sh $FM_MYAPP nogui"
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) runCommand "echo $'\033]30;' $FM_MYAPP Post Processing $'\007' "
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) runCommand "source $INSTALLATION_DIR/salome/salome_prerequisites.sh && source $INSTALLATION_DIR/salome/salome_modules.sh"
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) sendText " # Type \"paraview\" in order to visualize .med format solutions "
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
      echo "main= " $PWD
    CreateApp $1 dbg
     
    
  else
    # single tab
    CreateApp $1 $OPTION 
  fi
