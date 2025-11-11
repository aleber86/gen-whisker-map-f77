#set terminal pngcairo
#set output "time.png"
set size square
set xlabel "Condiciones iniciales"
set ylabel "Tiempo [s]"
set key outside center bottom
set grid lw 1.5
f(x) = a*x 
g(x) = d*x 

fit f(x) "time_data_total.time"  u 1:(($2)) via a 
fit g(x) "time_total.time" u 1:(($2)) via d

#set xrange[0:18]
plot "time_data_total.time" u 1:(($2)) pt 7 ps 2 title "GPU",\
"time_total_plus.time" u 1:(($2)) pt 7 ps 2 title "CPU",\
f(x) lw 2 lc rgb "blue" title "GPU",\
g(x) lw 2 lc rgb "red" title "CPU"
