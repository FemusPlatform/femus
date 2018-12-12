red=`tput setaf 9`; bold=`tput bold `; green=`tput setaf 10`; reset=`tput sgr0`; NC=`tput sgr0`
blue=`tput setaf 14`
# SET OF FUNCTIONS FOR USING FEMUS CODE

source functions.sh

function  command_exists () {
    type "$1" &> /dev/null ;
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
   select_tutorial_case $TUTORIAL_CLASS
      
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

