set xr [6520:6590]
set yr [0:9.5e-12]
set term png
set output 'fit_6881.png'
plot 'cngc6881_3300s.0054.txt' w l, 'fits.txt' w l
