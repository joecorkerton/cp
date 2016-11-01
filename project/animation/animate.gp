set term pngcairo enhanced font "sans, 40" size 3200,1800
set output 'h10.png'
set xlabel ''
set ylabel ''
set title 'Snapshot of Spins in Thermal Equilibrium when T = 3.0 J/k_B'
unset key
unset colorbox
unset ytics
unset xtics
set size square
set palette defined (-1 "black", 1 "red")
plot 'h10.txt' u 2:1:3 w image
