red=`tput setaf 9`; bold=`tput bold `; green=`tput setaf 10`; reset=`tput sgr0`; NC=`tput sgr0`
blue=`tput setaf 14`; 
# SET OF FUNCTIONS FOR USING FEMUS CODE

function  command_exists () {
    type "$1" &> /dev/null ;
}

function femus_guide {
   femus_show_application_functions
   femus_show_configure_functions
   femus_show_compiling_functions
   femus_show_tutorial_functions
}

function femus_show_configure_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for configuring femus applications${NC}"
  echo "  ${red}femus_application_configure <method>${NC}: function to be called within application folder."
  echo "     Application name, path and method (opt or dbg for compiling in optimized or debug mode) are set"
}


function femus_application_configure () {

   echo "-------------------------------------------------------"
   echo "Configuring application"
   
   unset APP_PATH
   unset METHOD
   unset FM_MYAPP
   
   export APP_PATH
   export METHOD
   export FM_MYAPP
   while :
   do
    case "$1" in
      opt | OPT) # optimized
          METHOD='opt'
          break
          ;;
      dbg | DBG) # debug
          METHOD='dbg'
          break
          ;;    
      *)  # unrecognized
          echo "Not available compilation method "$1
          echo "Rerun function femus_application_configure with valid method "
          echo "opt or dbg"
          return;
          break
          ;; 
    esac
   done 
   
   while :
   do
    case "$2" in
      "") # Application here
          APP_PATH=$PWD
          break
          ;;
      *)  # Application ad path $1
          APP_PATH=$2
          break
          ;; 
    esac
   done
   
   cd $APP_PATH
   cd ../
   local app_pre=$PWD
   cd $APP_PATH
   
   export path_len=${#app_pre}
   FM_MYAPP=${APP_PATH:$((path_len +1))}
   
   export HDF5_USE_FILE_LOCKING="FALSE"
   
   # Final print ===================================
   echo "Application path is " $APP_PATH 
   echo "Application name is "$FM_MYAPP
   echo "Method is "$METHOD
   echo "Application configured"
   # ==========================================================
   return;
}

# ==================================================================================
#                      COMPILE LIBRARY AND STANDARD APPLICATIONS
# ==================================================================================
function femus_show_compiling_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for compiling femus and  standard applications${NC}"
  echo "  ${red}femus_FEMuS_compile_lib${NC}: compile femus source files as separate libraries for 2D and 3D geometries"
  echo "  ${red}femus_gencase_compile_lib${NC}: compile standard gencase_2d and gencase_3d applications for QUAD 4/9 and HEX 8/27"
  echo "     2D and 3D meshes"
}



function femus_gencase_compile_lib {

  OLD_METHOD=$METHOD
  export ACTUAL_DIR=$PWD
  
  export MAIN_GENCASE_DIR=$FEMUS_DIR/applications/gencase/
  cd $MAIN_GENCASE_DIR
  
  # 2D gencase
  echo "${red}Now compiling standard gencase application for 2D geometries"
  echo "Standard gencase for quadrilateral elements${NC}"
  cd gencase_2d
  femus_application_configure opt 
  make clean
  make -j2
  if [ -f "$FEMUS_DIR/bin/gencase_2d" ]; then
    rm $FEMUS_DIR/bin/gencase_2d 
  fi
  ln -s $MAIN_GENCASE_DIR/gencase_2d/gencase_2d $FEMUS_DIR/bin/gencase_2d
  
  # 3D gencase
  echo "${red}Now compiling standard gencase application for 3D geometries"
  echo "Standard gencase for hexahedral elements${NC}"
  cd ../gencase_3d
  femus_application_configure opt 
  make clean
  make -j2
  if [ -f "$FEMUS_DIR/bin/gencase_3d" ]; then
    rm $FEMUS_DIR/bin/gencase_3d 
  fi
  ln -s $MAIN_GENCASE_DIR/gencase_3d/gencase_3d $FEMUS_DIR/bin/gencase_3d
  
  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
  
}

function femus_FEMuS_compile_lib {

  OLD_METHOD=$METHOD
  export ACTUAL_DIR=$PWD
  echo "Compiling femus library for 2D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus2D
  femus_application_configure opt
  # removing libfemus_2d.so
  make clean
  # removing object files from femus/src
  make src_clean
  make 
  
  echo "Compiling femus library for 3D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus3D
  femus_application_configure opt
  # removing libfemus_2d.so
  make clean
  # removing object files from femus/src
  make src_clean
  make 
  
  # cleaning after lib building
  make src_clean
  
  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
}

# ==================================================================================
#                              RUN APPLICATIONS
# ==================================================================================

function femus_show_application_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for running applications${NC}"
  echo "  ${red}femus_FEMuS_run <np>${NC}: compile actual application using np procs, and then run, using np procs"
  echo "  ${red}femus_gencase_run <np>${NC}: compile gencase using np procs, and then run, using np procs"
  echo "  ${red}femus_gencase_run_lib2D <np>${NC}: compile standard 2D gencase using np procs, and then run, using np procs"
  echo "  ${red}femus_gencase_run_lib3D <np>${NC}: compile standard 3D gencase using np procs, and then run, using np procs"
  echo "  ${red}femus_interpolator_run <np>${NC}: compile with np procs and then run, using np procs"
}
# =================================================================================
function femus_FEMuS_run {

   if [  "$1" == ""  ]; then pp="1" ;    else    pp=$1;    fi
   make -j$pp
   
   if [ "$?" != 0 ]; then
      echo "${red}${bold}=============================================="
      echo "             COMPILATION ERROR"
      echo "==============================================${reset}"
   else     
     mpiexec -np $pp $FM_MYAPP-$METHOD 2> messages.log
   fi
}

# =================================================================================
function femus_gencase_run {

   if [  "$1" == ""  ]; then pp="1" ;    else    pp=$1;    fi

   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     if [ ! -f "$FEMUS_DIR/bin/gencase" ]; then
       echo gencase not found make .....
       make gencase -j$pp
     fi
     echo run gencase .. $FEMUS_DIR/bin/gencase
     mpiexec -np $pp $FEMUS_DIR/bin/gencase 
     
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make gencase   "
     echo "      femus_gencase_run <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}
# =================================================================================
function femus_gencase_run_lib2D {
if [  "$1" == ""  ]; then pp="1" ;    else    pp=$1;    fi
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     
     if [ ! -f "$FEMUS_DIR/bin/gencase_2d" ]; then
       femus_gencase_compile_lib
     fi
     
     mpiexec -np $pp gencase_2d 2> messages_gencase.log
     
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make gencase   "
     echo "      femus_gencase_run <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}
# =================================================================================
function femus_gencase_run_lib3D {
if [  "$1" == ""  ]; then pp="1" ;    else    pp=$1;    fi
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     
     if [ ! -f "$FEMUS_DIR/bin/gencase_3d" ]; then
       femus_gencase_compile_lib
     fi
     
     mpiexec -np $pp gencase_3d 2> messages_gencase.log
     
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make gencase   "
     echo "      femus_gencase_run <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}
# =================================================================================
function femus_interpolator_run {
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     mpiexec -np $1 Interpolator 2> messages_interpolator.log
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make interpolator   "
     echo "      femus_interpolator_run <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}

# ==================================================================================
#                                  TUTORIALS
# ==================================================================================

function femus_show_tutorial_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for tutorial applications${NC}"
  echo "  ${red}femus_show_tutorials${NC}: show all available tutorials"
  echo "  ${red}femus_tutorial_run <path>${NC}: guided choice of tutorial for execution in <path> directory."
  echo "     If <path> is not given then tutorial is run in actual position"
}


function femus_show_tutorials () {
  if command_exists tree ; then
    echo "Tutorial folder structure:"
    echo "${red} Tutorial ${NC} - ${green} Class ${NC} - ${blue} Case ${NC}"
    tree -d -L 2 $FEMUS_DIR/tutorials
  else
    echo "${red} tree utility not present ${NC}"
    echo "Please install tree and re-run command"
    return;
  fi
  return;
}

function femus_tutorial_run () {
while :
do
    case "$1" in
      -h | --help)
          display_help  # Call your function
          exit 0
          ;;
      "")  # No more options
          echo "Tutorial will be run here"
          echo "Confirm your choice"
          select confirm in y n ;
          
          do
           case "$confirm" in
            y)
             break
             ;;
            n)
             return
             ;;
           esac
          done
          
          echo "__________________________________________"
          femus_tutorial_copy $PWD
          break
          ;;
      *)  # No more options
          echo "Tutorial will be run in the following path"
          echo "path: "$1
          echo "__________________________________________"
          femus_tutorial_copy $1
          break
          ;; 
    esac
done
return;
}




function femus_tutorial_run_all () {

export TUTORIAL_RUN
export TUTORIAL_HOME=$FEMUS_DIR/tutorials/

while :
do
    case "$1" in
      -h | --help)
          display_help  # Call your function
          exit 0
          ;;
      "")  # No more options
          TUTORIAL_RUN=$PWD
          echo "Tutorial will be run here"
          echo "Confirm your choice"
          select confirm in yes no ;
          
          do
           case "$confirm" in
            yes)
             break
             ;;
            no)
             return
             ;;
           esac
          done
          break
          ;;
      *)  # No more options
          echo "Tutorials will be run in the following path"
          echo "path: "$1
          TUTORIAL_RUN=$1
          break
          ;; 
    esac
done

TUTORIAL_RUN=$TUTORIAL_RUN/allTutorials
export TUTORIAL_LOG=$TUTORIAL_RUN/tutorialsLog.log

if [ -f "$TUTORIAL_LOG" ]; then
  rm $TUTORIAL_LOG
fi  
touch $TUTORIAL_LOG

export CLASSES=$(ls -l $TUTORIAL_HOME | grep ^d | awk '{print $9}')

if [ -d $TUTORIAL_RUN ];  then
  rm -r $TUTORIAL_RUN
fi  

mkdir $TUTORIAL_RUN
for class in $CLASSES; do

    echo "==========================================================================================">> $TUTORIAL_LOG
    echo "  NOW RUNNING $class CLASS TUTORIALS">> $TUTORIAL_LOG
    echo "==========================================================================================">> $TUTORIAL_LOG

    cd $TUTORIAL_RUN
    mkdir $TUTORIAL_RUN/$class
    unset CASES
    export CASES=$(ls -l $TUTORIAL_HOME/$class | grep ^d | awk '{print $9}')
    for tutorial in $CASES; do
        echo "******************************************************************************************">> $TUTORIAL_LOG
        echo "  NOW RUNNING $tutorial CASE ">> $TUTORIAL_LOG
        echo "******************************************************************************************">> $TUTORIAL_LOG
       unset tutorial_path
       export tutorial_path=$TUTORIAL_RUN/$class/$tutorial
       cp -r $TUTORIAL_HOME/$class/$tutorial $tutorial_path
       cd $tutorial_path
       femus_application_configure opt  >> $TUTORIAL_LOG
       source runTest.sh 
       cd $TUTORIAL_RUN/$class
       echo >> $TUTORIAL_LOG
    done    
    echo >> $TUTORIAL_LOG
    echo >> $TUTORIAL_LOG
done

cd $TUTORIAL_RUN

echo "ALL TUTORIALS: COMPLETED RUN"

}


function femus_tutorial_copy () {
   
   while :
do
    case "$1" in
      -h | --help)
          display_help  # Call your function
          exit 0
          ;;
      "")  # No more options
          TUTORIAL_CP=$PWD
          echo "Tutorial will be copied here"
          echo "Confirm your choice"
          select confirm in yes no ;
          
          do
           case "$confirm" in
            yes)
             break
             ;;
            no)
             return
             ;;
           esac
          done
          break
          ;;
      *)  # No more options
          echo "Tutorials will be copied in the following path"
          echo "path: "$1
          TUTORIAL_CP=$1
          break
          ;; 
    esac
done
   cd $TUTORIAL_CP
   echo "Available tutorials are"
   femus_show_tutorials
   
   export TUTORIAL_PATH=$FEMUS_DIR/tutorials/
         
   export TUTORIAL_CLASS
   export TUTORIAL_CASE
   femus_tutorial_class_select
#    femus_tutorial_select_case $TUTORIAL_CLASS
      
   echo "Selected test is $TUTORIAL_CASE of class $TUTORIAL_CLASS"
   
   cp -r $FEMUS_DIR/tutorials/$TUTORIAL_CLASS/$TUTORIAL_CASE $TUTORIAL_CP
   
   cd $TUTORIAL_CASE
   
   echo "-------------------------------------------------------"
   echo "Tutorial has been copied"
   echo "To run tutorial please execute"
   echo "femus_application_configure <Method> (opt or dbg)"
   echo "and then"
   echo "source runTest.sh"
   echo "-------------------------------------------------------"
   
   
   return;
}

function femus_tutorial_class_select () {

 export CLASSES=$(ls -l $FEMUS_DIR/tutorials | grep ^d | awk '{print $9}')
  
 echo "Choose tutorial class: " 
 select class in $CLASSES ;
  do
   case "$class" in
    *) echo $class " class chosen "
    TUTORIAL_CLASS=$class
    break;;
   esac
  done

  femus_tutorial_select_case $TUTORIAL_CLASS
}

function femus_tutorial_select_case () {
 export CASES=$(ls -l $FEMUS_DIR/tutorials/$1 | grep ^d | awk '{print $9}')
  
 echo "Choose tutorial case for class $TUTORIAL_CLASS: " 
 select tut_case in $CASES back_to_classes ;
  do
   case "$tut_case" in
    back_to_classes)
       femus_tutorial_class_select
    break;;   
    *) echo $tut_case " case chosen"
    TUTORIAL_CASE=$tut_case
    break ;;
   esac
  done

}
