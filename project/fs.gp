set term pngcairo enhanced font "sans, 40" size 3200,1800
set output 'fs.png'
set xlabel '1 / (N * N)^2'
set ylabel 'Critical Temperature [J/k_B]'
set title 'Plot of Critical Temperature against (Lattice Size)^{-2} '
unset key
set yrange [2.25:2.55]
set xrange [0:0.045]
f(x) = 4.2*x + 2.332
plot 'fs.txt' u 2:4 w p ps 10, f(x)
#f(x) = a*x + b
#fit f(x) 'fs.txt' u 2:4  via a, b
