set title 'Monte-Carlo Integration'
set xlabel 'log (Number of sample points)'
set ylabel 'log (Error)'
set logscale xy
set nokey
#y(x) = a*x**(-0.5)
#a = 150.0
#fit y(x) 'output.txt' via a
plot 'output.txt' u 1:3, 35*x**(-0.5) 
set term pngcairo size 3200,1800
set output 'monte-carlo.png'
replot
pause -1
