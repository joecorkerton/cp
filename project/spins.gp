set term pngcairo enhanced font "sans, 40" size 3200,1800
set output 'H.png'
set xlabel ''
set ylabel ''
set title 'Snapshot of Spins in Thermal Equilibrium when T = 2.27 J/k_B with an External Field H = 1 T'
unset key
unset colorbox
unset ytics
unset xtics
set size square
set palette defined (-1 "black", 1 "red")
plot 'H.txt' u 2:1:3 w image
