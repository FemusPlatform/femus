
# ====================================================================
# Function check with 3 arguments: 
# $1= check_status
# $2=case
# $3=file to write (test.log) label 
function CheckStatus {
  if [ $1 == 0 ]; then
  # if  check_status=0
    echo "   OK: case= " $2  >> $3
    echo "   OK: case= " $2 
    echo " ========================================"
  else
  # if  check_status=0
    echo "   Wrong exit: case= " $2 " with status " $1 >> $3  
    echo " ========================================"
    echo "   Fail to run: case= " $2 " with status " $1  
     echo " ========================================"
     exit $1
  fi
}

# ====================================================================
# Function check with 3 arguments: 
# $1= simul -> label RESU
# $2=case
# $3= processors number

function SetTest {

SIMU=$1
CASE=$2

echo "---------------------------------------------------------------------">> test.log
echo " Now executing test " $SIMU >> test.log
echo "---------------------------------------------------------------------">> test.log
echo " Now executing test" $SIMU 
echo "---------------------------------------------------------------------">> test.log
echo "---------------------------------------------------------------------">> test.log
export exit_status=1

if [ -d "DATA" ]; then
  rm -fr DATA
fi
if [ -d "SRC" ]; then
  rm -fr SRC
fi
if [ -d "MESH" ]; then
  rm -fr MESH
fi
  mkdir DATA
  mkdir SRC
  mkdir MESH
# echo " Copy files from" $SIMU 
cp  test/$SIMU/DATA/*.*      ./DATA/
cp  test/$SIMU/SRC/*.*       ./SRC/
cp  test/$SIMU/MESH/*.*      ./MESH/

 
 
# Note:  $? it's the exit status of the last command
make resu_clean
if [ -d "RESU_MED" ]; then
  rm -fr RESU_MED
fi
make -j3 gencase > test.log
exit_status=$?
CheckStatus $exit_status " make gencase in "$1 test.log


mpiexec -np $3 ../gencase/gencase-opt > test.log
exit_status=$?
CheckStatus $exit_status " execution gencase in "$1 test.log


make -j3 > test.log
exit_status=$?
CheckStatus $exit_status " make femus_test-opt in "$1 test.log


mpiexec -np $3 femus_test-opt 
exit_status=$?
# CheckStatus $exit_status " execution NS_test-opt in "$1 test.log
echo " "
echo " ========================================"
echo " End test" $SIMU 


if [ -d "RESU_"$SIMU ]; then
  rm -fr RESU_$SIMU
fi


mkdir RESU_$SIMU
mkdir RESU_$SIMU/RESU
mkdir RESU_$SIMU/RESU_MED
cp ./RESU/*.* ./RESU_$SIMU/RESU
cp ./RESU_MED/*.* ./RESU_$SIMU/RESU_MED

make resu_clean
make src_clean

# echo >> test.log
# echo >> test.log

}

# ====================================================
#                     main
# ====================================================
# comments -> test.log
if [ -f "test.log" ]; then
  rm test.log
fi
echo " ========================================"
echo " ========================================"



# ===================================================
#  base
# ===================================================

 SetTest     test_nse   test_nse  1
 SetTest     test_ctrlfsi  test_ctrlfsi  1

echo " ========================================"
