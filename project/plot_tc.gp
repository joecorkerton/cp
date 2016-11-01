set term pngcairo enhanced font "sans, 40" size 3200,1800
set output 'tc_hc.png'
#set xlabel 'Temperature [J/k_B]'
#set ylabel '|Magnetisation| [{/Symbol m}]'
#set title 'Plot of Average Magnetisation Per Site for a Range of Different Initial Temperatures'
#set xlabel 'Temperature [J/k_B]'
#set ylabel 'Energy [J]'
#set title 'Plot of Average Energy Per Site for a Range of Different Initial Temperatures'
set xlabel 'Temperature [J/k_B]'
set ylabel 'Heat Capacity [k_B^2 K/J]'
set title 'Plot of Average Heat Capacity Per Site for a Range of Different Initial Temperatures'
unset key
plot 'tc.txt' u 2:5 w linesp 

