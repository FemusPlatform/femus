  red=`tput setaf 1`
  bold=`tput bold `
  green=`tput setaf 2`
  reset=`tput sgr0`

# copy standard Makefile for applications
cp $FEMUS_DIR/solvers/Makefile .

# create directories for results
mkdir RESU
mkdir RESU_MED

# copy standard solvers into local SRC folder
cp $FEMUS_DIR/solvers/MGSolverNS/*.h SRC/
cp $FEMUS_DIR/solvers/MGSolverNS/*.C SRC/

echo "Compiling application "
if [ -f $FEMUS_DIR/lib/libfemus_3d.so ]; then
  echo "FEMuS 3D library found "
  echo "Compiling with library "
  make withlib_3d -j2
else
  echo "FEMuS 3D library NOT found "
  echo "Compiling without library "
  make -j2
fi  
  
if [ "$?" != 0 ]; then
      echo "${red}${bold}=============================================="
      echo "             APPLICATION COMPILATION ERROR"
      echo "==============================================${reset}"
      return
fi     

echo "Compiling gencase "
make gencase

if [ "$?" != 0 ]; then
      echo "${red}${bold}=============================================="
      echo "             GENCASE COMPILATION ERROR"
      echo "==============================================${reset}"
      return
fi    

runGencase 1

mpiexec -np 1 $FM_MYAPP-opt 2> messages.log

