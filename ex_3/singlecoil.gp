set term pngcairo size 3200,1800
set output 'Axial_field_diff.png'
set xlabel 'x /m'
set ylabel 'Magnetic Field /T'
set title 'Plot of difference between calculated B field and analytic result against displacement along the x axis for a single coil'
unset key
set format y "%.1s*10^{%S}"
plot 'singlecoilx.txt' u 1:4 w linesp

