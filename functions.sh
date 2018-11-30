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
