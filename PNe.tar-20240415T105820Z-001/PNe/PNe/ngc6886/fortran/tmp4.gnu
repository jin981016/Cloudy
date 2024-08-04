set xr [6520:6600]
#set xr [6550:6580]
set yr [0:1e-11]
#set yr [0:2e-14]
set term png
set output 'fits_6886.png'
plot 'ngc6886_2400s.0051.txt' w l, 'fits.txt' w l
