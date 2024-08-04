reset
set encoding iso_8859_1
set term png size 800,600
set output 'alpha_6741.png'

set xr[6520:6590]
#set xr[6540:6550]
#set yr[0:6.5e-11]
set yr[0:1.5e-13]
set xlabel 'Wavelength({\305})' font 'arial,13'
set ylabel 'Flux(erg s^{-1} {\305}^{-1})' font 'arial,13'
set label 'Raman He II' at 6538, 7.1e-14 font 'arial,15' front
set label '6545' at 6541, 6.5e-14 font 'arial,15' front
#set label 'H{/Symbol a}' at 6567,1.2e-13 font 'arial,15'
set label '[N II]' at 6530, 1.1e-13 font 'arial,15'
set arrow 1 from 6543, 6e-14 to 6545,4.8e-14
set title 'UVES spectrum of NGC6741' font 'arial,15'
plot 'cngc6741_3600s.0055.txt' w l lw 0.5 lt 1 lc "black" notitle,'fit_a.txt' w l lw 2.5 lt 3 lc "red" notitle

