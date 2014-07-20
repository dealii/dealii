set data style lines
set hidden3d

set xrange [-10:10]
set yrange [-2.2:2.2]
set zrange [-1:8]
set xlabel "x"
if ($0==1) set ylabel "u"
if ($0==2) set ylabel "y"
if ($0==2) set zlabel "u"

set terminal png
set output "solution-$0d-$1.png" 

if ($0==1) pl "solution-$0d-$1.gpl"
if ($0==2) spl "solution-$0d-$1.gpl"

#set terminal x11
#replot
#pause 5
