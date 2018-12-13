red=`tput setaf 9`; bold=`tput bold `; green=`tput setaf 10`; reset=`tput sgr0`; NC=`tput sgr0`
blue=`tput setaf 14`; 
# SET OF FUNCTIONS FOR USING FEMUS CODE

function  command_exists () {
    type "$1" &> /dev/null ;
}

function femusGuide {
   show_application_functions
   show_configure_functions
   show_compiling_functions
   show_tutorial_functions
}

function show_configure_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for configuring femus applications${NC}"
  echo "  ${red}configureApplication <method>${NC}: function to be called within application folder."
  echo "     Application name, path and method (opt or dbg for compiling in optimized or debug mode) are set"
}

function configureApplication () {
   echo
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
          echo "Rerun function configure_application with valid method "
          echo "opt or dbg"
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
   
   echo "Application path is " $APP_PATH 
   echo "Application name is "$FM_MYAPP
   echo "Method is "$METHOD
   echo "-------------------------------------------------------"
   echo "Application configured"
   echo "-------------------------------------------------------"
   
   return;
}

# ==================================================================================
#                      COMPILE LIBRARY AND STANDARD APPLICATIONS
# ==================================================================================
function show_compiling_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for compiling femus and  standard applications${NC}"
  echo "  ${red}compileLibrary${NC}: compile femus source files as separate libraries for 2D and 3D geometries"
  echo "  ${red}compileGencase${NC}: compile standard gencase_2d and gencase_3d applications for QUAD 4/9 and HEX 8/27"
  echo "     2D and 3D meshes"
}



function compileGencase {
  export ACTUAL_DIR=$PWD
  
  export MAIN_GENCASE_DIR=$FEMUS_DIR/applications/gencase/
  cd $MAIN_GENCASE_DIR
  
  # 2D gencase
  echo "${red}Now compiling standard gencase application for 2D geometries"
  echo "Standard gencase for quadrilateral elements${NC}"
  cd gencase_2d
  configureApplication opt 
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
  configureApplication opt 
  make clean
  make -j2
  if [ -f "$FEMUS_DIR/bin/gencase_3d" ]; then
    rm $FEMUS_DIR/bin/gencase_3d 
  fi
  ln -s $MAIN_GENCASE_DIR/gencase_3d/gencase_3d $FEMUS_DIR/bin/gencase_3d
  
  cd $ACTUAL_DIR

}

function compileLibrary {

  echo "Compiling femus library for 2D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus2D
  configure_application opt
  # removing libfemus_2d.so
  make clean
  # removing object files from femus/src
  make src_clean
  make 
  
  echo "Compiling femus library for 3D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus3D
  configure_application opt
  # removing libfemus_2d.so
  make clean
  # removing object files from femus/src
  make src_clean
  make 
  
  # cleaning after lib building
  make src_clean

}

# ==================================================================================
#                              RUN APPLICATIONS
# ==================================================================================

function show_application_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for running applications${NC}"
  echo "  ${red}runFEMuS <np>${NC}: compile actual application using np procs, and then run, using np procs"
  echo "  ${red}runGencase <np>${NC}: compile gencase using np procs, and then run, using np procs"
  echo "  ${red}runGencase2D <np>${NC}: compile standard 2D gencase using np procs, and then run, using np procs"
  echo "  ${red}runGencase3D <np>${NC}: compile standard 3D gencase using np procs, and then run, using np procs"
  echo "  ${red}runInterpolator <np>${NC}: compile with np procs and then run, using np procs"
}

function runFEMuS {
   make -j$1
   
   if [ "$?" != 0 ]; then
      echo "${red}${bold}=============================================="
      echo "             COMPILATION ERROR"
      echo "==============================================${reset}"
   else     
     mpiexec -np $1 $FM_MYAPP-opt 2> messages.log
   fi
}

function runGencase {
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     if [ ! -f "$FEMUS_DIR/bin/gencase" ]; then
       make gencase
     fi
     mpiexec -np $1 gencase 2> messages_gencase.log
     
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make gencase   "
     echo "      runGencase <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}

function runGencase2D {
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     
     if [ ! -f "$FEMUS_DIR/bin/gencase_2d" ]; then
       compileGencase
     fi
     
     mpiexec -np $1 gencase_2d 2> messages_gencase.log
     
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make gencase   "
     echo "      runGencase <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}

function runGencase3D {
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     
     if [ ! -f "$FEMUS_DIR/bin/gencase_3d" ]; then
       compileGencase
     fi
     
     mpiexec -np $1 gencase_3d 2> messages_gencase.log
     
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make gencase   "
     echo "      runGencase <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}

function runInterpolator {
   if [ -d  $FEMUS_DIR/bin ]; then
     echo "found bin directory "
     mpiexec -np $1 Interpolator 2> messages_interpolator.log
   else
     echo "${red} ======================================================= "
     echo " UNABLE TO FIND /USER_APPL/bin DIRECTORY  "
     echo " CREATE bin DIRECTORY, THEN TYPE "
     echo "      make interpolator   "
     echo "      runInterpolator <n proc> "
     echo " ======================================================= ${reset}"
   fi   
}

# ==================================================================================
#                                  TUTORIALS
# ==================================================================================

function show_tutorial_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for tutorial applications${NC}"
  echo "  ${red}showTutorials${NC}: show all available tutorials"
  echo "  ${red}runTutorial <path>${NC}: guided choice of tutorial for execution in <path> directory."
  echo "     If <path> is not given then tutorial is run in actual position"
}


function showTutorials () {
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

function runTutorial () {
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
          copy_tutorial $PWD
          break
          ;;
      *)  # No more options
          echo "Tutorial will be run in the following path"
          echo "path: "$1
          echo "__________________________________________"
          copy_tutorial $1
          break
          ;; 
    esac
done
return;
}

function copy_tutorial () {
   echo "Available tutorials are"
   cd $1
   showTutorials
   
   export TUTORIAL_PATH=$FEMUS_DIR/tutorials/
         
   export TUTORIAL_CLASS
   export TUTORIAL_CASE
   select_tutorial_class
#    select_tutorial_case $TUTORIAL_CLASS
      
   echo "Selected test is $TUTORIAL_CASE of class $TUTORIAL_CLASS"
   
   cp -r $FEMUS_DIR/tutorials/$TUTORIAL_CLASS/$TUTORIAL_CASE $1
   
   cd $TUTORIAL_CASE
   
   echo "-------------------------------------------------------"
   echo "Tutorial has been copied"
   echo "To run tutorial please execute"
   echo "configureApplication <Method> (opt or dbg)"
   echo "and then"
   echo "source runTest.sh"
   echo "-------------------------------------------------------"
   
   
   return;
}

function select_tutorial_class () {

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

  select_tutorial_case $TUTORIAL_CLASS
}

function select_tutorial_case () {
 export CASES=$(ls -l $FEMUS_DIR/tutorials/$1 | grep ^d | awk '{print $9}')
  
 echo "Choose tutorial case for class $TUTORIAL_CLASS: " 
 select tut_case in $CASES back_to_classes ;
  do
   case "$tut_case" in
    back_to_classes)
       select_tutorial_class
    break;;   
    *) echo $tut_case " case chosen"
    TUTORIAL_CASE=$tut_case
    break ;;
   esac
  done

}
