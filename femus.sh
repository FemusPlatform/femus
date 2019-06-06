red=`tput setaf 9`; 
bold=`tput bold `; 
green=`tput setaf 10`; 
reset=`tput sgr0`;
NC=`tput sgr0`;
blue=`tput setaf 14`; 
# SET OF FUNCTIONS FOR USING FEMUS CODE

function  command_exists () {
    type "$1" &> /dev/null ;
}

function femus_guide {
   intern_femus_show_application_functions
   intern_femus_show_configure_functions
   intern_femus_show_compiling_functions
   intern_femus_show_tutorial_functions
   intern_femus_repo_functions
}
git commit -m "intern prefix added to some femus functions. femus + <tab> now shows only user available functions "
function intern_femus_repo_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for using femus and github repo${NC}"
  echo "  ${red}femus_install_repo_format${NC}: the function creates an executable pre-commit"
  echo "     hook file inside femus/.git/ folder. The command is used to automatically format"
  echo "     source files, added with \"git add\" function and before \"git commit\", using clang"
}

function intern_femus_show_configure_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for configuring femus applications${NC}"
  echo "  ${red}femus_application_configure <method>${NC}: function to be called within application folder."
  echo "     Application name, path and method (opt or dbg for compiling in optimized or debug mode) are set"
  echo "  ${red}femus_link_solver_files${NC}: function for linking solver files inside femus/solvers/ln_solvers folder"
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
function intern_femus_show_compiling_functions {
  echo "------------------------------------------------------------"
  echo "${green}List of functions for compiling femus and  standard applications${NC}"
  echo "  ${red}femus_FEMuS_compile_lib${NC}: compile femus source files as separate libraries for 2D and 3D geometries"
  echo "  ${red}femus_gencase_compile_lib${NC}: compile standard gencase_2d and gencase_3d applications for QUAD 4/9 and HEX 8/27"
  echo "     2D and 3D meshes"
}



function femus_gencase_compile_lib {

  if [  "$1" == ""  ]; then np="1" ;    else    np=$1;    fi
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
  make -j$np
  
  # 3D gencase
  echo "${red}Now compiling standard gencase application for 3D geometries"
  echo "Standard gencase for hexahedral elements${NC}"
  cd ../gencase_3d
  femus_application_configure opt 
  make clean
  make -j$np

  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
  
}

function femus_turbulence_compile_lib_opt {

  if [  "$1" == ""  ]; then np="1" ;    else    np=$1;    fi
  OLD_METHOD=$METHOD
  export ACTUAL_DIR=$PWD
  echo "Compiling femus turb library for 2D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_turb2D
  femus_application_configure opt
  # removing libturb_2d.so
  make clean
  make -j$np
  
  echo "Compiling femus turb library for 3D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_turb3D
  femus_application_configure opt
  # removing libturb_3d.so
  make clean
  make -j$np
  
  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
}

function femus_turbulence_compile_lib_dbg {

  if [  "$1" == ""  ]; then np="1" ;    else    np=$1;    fi
  OLD_METHOD=$METHOD
  export ACTUAL_DIR=$PWD
  echo "Compiling femus turb library for 2D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_turb2D
  femus_application_configure dbg
  # removing libturb_2d.so
  make clean
  # removing object files from femus/src
  make -j$np
  
  echo "Compiling femus turb library for 3D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_turb3D
  femus_application_configure dbg
  # removing libturb_3d.so
  make clean
  # removing object files from femus/src
  make -j$np
  
  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
}

function femus_FEMuS_compile_lib_opt {

  if [  "$1" == ""  ]; then np="1" ;    else    np=$1;    fi
  OLD_METHOD=$METHOD
  export ACTUAL_DIR=$PWD
  echo "Compiling femus library for 2D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus2D
  femus_application_configure opt
  # removing libfemus_2d.so
  make clean
  make -j$np
  
  echo "Compiling femus library for 3D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus3D
  femus_application_configure opt
  # removing libfemus_2d.so
  make clean
  make -j$np
  
  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
}

function femus_FEMuS_compile_lib_dbg {

  if [  "$1" == ""  ]; then np="1" ;    else    np=$1;    fi
  OLD_METHOD=$METHOD
  export ACTUAL_DIR=$PWD
  echo "Compiling femus library for 2D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus2D
  femus_application_configure dbg
  # removing libfemus_2d.so
  make clean
  # removing object files from femus/src
  make src_clean
  make -j$np
  
  echo "Compiling femus library for 3D geometry "
  cd $PLAT_CODES_DIR/femus/applications/lib_femus3D
  femus_application_configure dbg
  # removing libfemus_2d.so
  make clean
  # removing object files from femus/src
  make src_clean
  make -j$np
  
  # cleaning after lib building
  make src_clean
  
  cd $ACTUAL_DIR
  export METHOD=$OLD_METHOD
  femus_application_configure $METHOD
}

function femus_compile_all {

  if [  "$1" == ""  ]; then np="1" ;    else    np=$1;    fi
  
  femus_FEMuS_compile_lib_dbg $np
  femus_FEMuS_compile_lib_opt $np
  femus_turbulence_compile_lib_dbg $np
  femus_turbulence_compile_lib_opt $np
  femus_gencase_compile_lib $np
}

# ==================================================================================
#                              RUN APPLICATIONS
# ==================================================================================

function intern_femus_show_application_functions {
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
     echo run gencase .. $FEMUS_DIR/bin/gencase-opt
     mpiexec -np $pp $FEMUS_DIR/bin/gencase-opt 
     
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

function intern_femus_show_tutorial_functions {
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
          intern_femus_tutorial_copy $PWD
          break
          ;;
      *)  # No more options
          echo "Tutorial will be run in the following path"
          echo "path: "$1
          echo "__________________________________________"
          intern_femus_tutorial_copy $1
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

    EXEC_TEST=1
    if [ "$class" == "DD_coupling" ]; then
       if [ ! -d $PLAT_CODES_DIR/dragondonjon ]; then
          EXEC_TEST=0
       fi  
    fi
    
    if [ "$EXEC_TEST" == "1" ]; then
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
    fi
done

cd $TUTORIAL_RUN

echo "ALL TUTORIALS: COMPLETED RUN"

}


function intern_femus_tutorial_copy () {

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
   
   # here class and case are chosen from user
   intern_femus_tutorial_class_select
      
   if [ "$TUTORIAL_CASE" == "all" ]; then
   
   TUTORIAL_RUN=$TUTORIAL_CP/$TUTORIAL_CLASS"_tut"
   export TUTORIAL_LOG=$TUTORIAL_RUN/tutorialsLog.log
   
   
   if [ -d $TUTORIAL_RUN ]; then
     rm -r $TUTORIAL_RUN
   fi     
     mkdir $TUTORIAL_RUN
     cd $TUTORIAL_RUN
     
     touch $TUTORIAL_LOG
     
       unset CASES
       export CASES=$(ls -l $TUTORIAL_HOME/$class | grep ^d | awk '{print $9}')
       for tutorial in $CASES; do
           echo "******************************************************************************************">> $TUTORIAL_LOG
           echo "  NOW RUNNING $tutorial CASE ">> $TUTORIAL_LOG
           echo "******************************************************************************************">> $TUTORIAL_LOG
          unset tutorial_path
          export tutorial_path=$TUTORIAL_RUN/$tutorial
          cp -r $TUTORIAL_HOME/$class/$tutorial $tutorial_path
          cd $tutorial_path
          femus_application_configure opt  >> $TUTORIAL_LOG
          source runTest.sh 
          cd $TUTORIAL_RUN
          echo >> $TUTORIAL_LOG
       done    
       echo >> $TUTORIAL_LOG
       echo >> $TUTORIAL_LOG
   
   else
      
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
   
   fi
   
   return;
}

function intern_femus_tutorial_class_select () {

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

  intern_femus_tutorial_select_case $TUTORIAL_CLASS
}

function intern_femus_tutorial_select_case () {
 export CASES=$(ls -l $FEMUS_DIR/tutorials/$1 | grep ^d | awk '{print $9}')
  
 echo "Choose tutorial case for class $TUTORIAL_CLASS: " 
 select tut_case in $CASES all back_to_classes ;
  do
   case "$tut_case" in
    back_to_classes)
       intern_femus_tutorial_class_select
    break;;   
    "all")
    echo " all tutorials of class $TUTORIAL_CLASS will be run"
    TUTORIAL_CASE=$tut_case
    break ;;
    *) echo $tut_case " case chosen"
    TUTORIAL_CASE=$tut_case
    break ;;
   esac
  done

}

function femus_link_solver_files () {

  unset SOLVERS
  export SOLVERS=$(ls  $PLAT_CODES_DIR/femus/solvers/ | grep "MG" )
  
  if [ -d $PLAT_CODES_DIR/femus/solvers/ln_solvers/ ]; then
    rm -r $PLAT_CODES_DIR/femus/solvers/ln_solvers/
  fi
  mkdir $PLAT_CODES_DIR/femus/solvers/ln_solvers/
  
  for solver in $SOLVERS; do
      unset SRC_FILES_TO_LINK
      unset HEADER_FILES_TO_LINK
      export SRC_FILES_TO_LINK=$(find $PLAT_CODES_DIR/femus/solvers/$solver/ -maxdepth 3 -name "*.C" )
      export HEADER_FILES_TO_LINK=$(find $PLAT_CODES_DIR/femus/solvers/$solver/ -maxdepth 3 -name "*.h" )
      for source_file in $SRC_FILES_TO_LINK; do
          ln -s $source_file $PLAT_CODES_DIR/femus/solvers/ln_solvers/
      done
      for header_file in $HEADER_FILES_TO_LINK; do
          ln -s $header_file $PLAT_CODES_DIR/femus/solvers/ln_solvers/
      done
  done 
  
}

function femus_install_repo_format (){

PRE_HOOK_FORMAT=$FEMUS_DIR/.git/hooks/pre-commit

if [ -f $PRE_HOOK_FORMAT ]; then
  rm $PRE_HOOK_FORMAT
fi

touch $PRE_HOOK_FORMAT
chmod 777 $PRE_HOOK_FORMAT

echo "#!/bin/bash

STYLE=\$(git config --get hooks.clangformat.style)
if [ -n \"\${STYLE}\" ] ; then
  STYLEARG=\"-style=\${STYLE}\"
else
  STYLEARG=\"\"
fi

format_file() {
  file=\"\${1}\"
  clang-format -i \${STYLEARG} \${1}
  git add \${1}
}

case \"\${1}\" in
  --about )
    echo \"Runs clang-format on source files\"
    ;;
  * )
    for file in \`git diff-index --cached --name-only HEAD\` ; do

    filename=\`basename \$file\`
    extention=\`sed 's/^\w\+.//' <<< \"\$filename\"\`

      case \"\$extention\" in
      \"C\" ) 
         format_file \"\${file}\"
         echo \"File \" \"\$file\" \" has been formatted using clang-format \"
         break ;;
      \"h\" ) 
         format_file \"\${file}\"
         echo \"File \" \"\$file\" \" has been formatted using clang-format \"
         break ;;
      * ) 
         echo \"File \" \"\$file\" \" added without changing anything \"
         git add \$file
         break ;;
      esac
    done
    ;;
esac ">> $PRE_HOOK_FORMAT

echo "Pre-commit hook file has been created"

}


function femus_residual(){

  tee >(grep -oP --line-buffered 'Cumulative Residual '$1' \K.*'  > residual_$1.dat)

}
