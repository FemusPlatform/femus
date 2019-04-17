set terminal wxt size 500,400


set style line 1 lt rgb "red" lw 2 
set style line 12 lt rgb "red" lw 1 

set style line 2 lt rgb "blue" lw 2
set style line 22 lt rgb "blue" lw 1


set xtics
set ytics


set title "Residual plot"
set key left top
plot "residual_lk.dat"  u 1 w l ls 1 title "lk 1", \
     "residual_lw.dat"  u 1 w l ls 2 title "lw 1", 
#      "residual_lk.dat"  u 2 w l ls 12 title "lk 0", \
#      "residual_lw.dat"  u 2 w l ls 22 title "lw 0", 
 
pause 10
reread

