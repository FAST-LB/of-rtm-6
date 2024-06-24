set multiplot layout 1, 1 title "Pressures" font ",14"

set tmargin 2

set xlabel "Time"
set ylabel "Pressure [Pa]"
unset key

plot 'postProcessing/p_1/0/p' title 'p_1' with lines, 'postProcessing/p_2/0/p' title 'p_2' with lines, 'postProcessing/p_3/0/p' title 'p_3' with lines, 'postProcessing/p_4/0/p' title 'p_4' with lines, 'postProcessing/p_5/0/p' title 'p_5' with lines, 'postProcessing/p_6/0/p' title 'p_6' with lines,

pause 2
reread 