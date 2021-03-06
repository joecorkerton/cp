set xlabel 'Position x /m'
set ylabel 'Relative amplitude/phase'
set title 'Diffraction pattern for D = 100 cm'
plot 'Diffraction_100.txt' u 1:2 title 'Amplitude' w lines , 'Diffraction_100.txt' u 1:3 title 'Phase' w lines
set term pngcairo size 3200,1800
set output 'slit_100.png'
replot
pause -1
