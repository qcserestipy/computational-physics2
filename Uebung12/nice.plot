set key off
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"

set terminal pdf
set output "out.pdf"
plot "autocorr.out" using 1:2 w lines lt 3 lc 7 ,"autocorr.out" using 1:3 w lines lt 1 lc 8,"autocorr.out" using 1:4 w lines lt 1 lc 9
plot "power_spec_A.out" using 1:2 w lines lt 1 lc 7
plot "power_spec_B.out" using 1:2 w lines lt 1 lc 7
plot "power_spec_C.out" using 1:2 w lines lt 1 lc 7
plot "power_spec_D.out" using 1:2 w lines lt 1 lc 7
