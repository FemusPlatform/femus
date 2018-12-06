# ==========================================================================================
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
########## $PLAT_DIR, given when installing femus; #############
########## otherwise, write it on your own                 #############
########################################################################

# ===========================================================================================================
  red=`tput setaf 1`
  bold=`tput bold `
  green=`tput setaf 2`
  reset=`tput sgr0`
# ===========================================================================================================


# function to run FEMuS 
# the function takes as input the number of processors
# a log file is printed by redirecting the output of std::cerr
source functions.sh


  
# ===========================================================================================================
# function CheckExistence: check first argument $1 as name application
# ===========================================================================================================
function CheckExistence {
   cd $CurDir
   if [ -d  USER_APPL/$1 ]; then
     # USER_APPL/$FM_MYAPP  exists
     echo "USER_APPL/"$1 " directory exists  " 
   else
     # USER_APPL/$FM_MYAPP does not exist
     echo "USER_APPL/"$1 " directory does not  exist  "
     echo "Do you want to create the USER_APPL/"$1  "application ? [yes/no]"
     read test_app
     if test "$test_app" = "yes"; then 
       cp -r template_appl/$APP_NAME/   $CurDir/USER_APPL/$1
     else
       return
     fi
   fi
     cd $CurDir/USER_APPL/$1
     echo "APPLICATION=" $1 "APP_NAME=" $APP_NAME "APP_VER=" $APP_VER "METHOD=" $METHOD "in USER_APPL/"$FM_MYAPP
     echo  "path " $PWD 
     export APP_PATH=$PWD 
}

# ===========================================================================================================
# function CreateApp: create application with name $1 with mode dbg,opt, etc....
# ===========================================================================================================
function CreateApp {
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
   cd $FEMUS_DIR/USER_APPL/gencase
   return;
  fi
  if test "$OPTION" = "interpolator"; then
   $PWD
   cd $FEMUS_DIR/USER_APPL/Interpolator
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
cd $FEMUS_DIR/template_appl 
ls
if [ -d  $APP_NAME ]; then
  echo $APP_NAME " is an application template - ok  "
  cd $FEMUS_DIR
else
  echo  " "$1" name is invalid: not VALID application template or _n is missing "
  cd $FEMUS_DIR
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
  export CurDir=$PWD
  echo "Current Directory " $CurDir
  export FEMUS_PATH=$PLAT_CODES_DIR/femus

  export LL=`expr length $PLAT_DIR`
  export LLL=`expr length $FEMUS_PATH`
  export L1=`expr $LLL - $LL`
  
  export FEMUS_NAME=`expr substr $FEMUS_PATH $(($LL +1)) $L1`

  echo "---------------------------------------------------------------"
  echo "FEMUS PATH   "      $FEMUS_PATH
  echo "INSTALLATION DIR   "$PLAT_DIR
  echo "FEMUS NAME   "      $FEMUS_NAME
  echo "---------------------------------------------------------------"
  
  echo "INSTALLAION_PATH=" $PLAT_DIR
  # environment MPI PETSC libmesh HDF5 MED Medcoupling ......

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
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var1)) runCommand "cd $CurDir && source $PLAT_DIR/plat_conf.sh && source $FEMUS_DIR/configure_femus.sh $FM_MYAPP nogui"
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
    var2=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession) 
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var2)) runCommand "cd $CurDir  && source $PLAT_DIR/plat_conf.sh && source $FEMUS_DIR/configure_femus.sh $FM_MYAPP gencase"
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
    var3=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession)   
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var3)) runCommand "cd $CurDir  && source $PLAT_DIR/plat_conf.sh && source $FEMUS_DIR/configure_femus.sh $FM_MYAPP interpolator"
    var4=$(qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW newSession) 
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) runCommand "cd $CurDir  && source $PLAT_DIR/plat_conf.sh && source $FEMUS_DIR/configure_femus.sh $FM_MYAPP nogui"
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) runCommand "echo $'\033]30;' $FM_MYAPP Post Processing $'\007' "
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) runCommand "source $PLAT_DIR/plat_conf.sh && source $PLAT_THIRD_PARTY_DIR/salome/salome_prerequisites.sh && source $PLAT_THIRD_PARTY_DIR/salome/salome_modules.sh"
    qdbus $KONSOLE_DBUS_SERVICE /Sessions/$(($var4)) sendText " # Type \"paraview\" in order to visualize .med format solutions "
    qdbus $KONSOLE_DBUS_SERVICE $KONSOLE_DBUS_WINDOW prevSession
      echo "main= " $PWD
    CreateApp $1 dbg
     
    
  else
    # single tab
    CreateApp $1 $OPTION 
  fi
