 if [ ! -d RESU/ ]; then
   mkdir RESU
 fi
 if [ ! -d RESU_MED/ ]; then
   mkdir RESU_MED
 fi
 
 make withlib_3d -j2
 
 mpiexec -np 1 gencase_3d
 mpiexec -np 1 neut-opt
