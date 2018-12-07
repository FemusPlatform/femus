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
make -j2

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

runFEMuS 1
