set term pngcairo size 3200,1800
set output 'small_angle.png'
set xlabel 'time /s'
set ylabel 'Amplitude'
set title 'Plot of displacement/ velocity over 1000 periods of oscillation'
set xrange [5800:6000]
plot 'output.txt' u 1:2 title 'theta' w lines, 'output.txt' u 1:3 title 'd(theta)/dt' w lines

