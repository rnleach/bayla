set datafile separator ","

fn = "simple_model.csv"

set palette rgb 33,13,10
set grid
set tics nomirror

set origin 0,0
set size 1,1
set multiplot layout 2,2 rowsfirst scale 0.9, 0.9

set title "Q"
plot fn using 1:2:4 with points palette pt 7 notitle
set title "P"
plot fn using 1:2:3 with points palette pt 7 notitle
set title "Credible Interval"
plot fn using 1:2:6 with points palette pt 7 notitle
set title "W"
plot fn using 1:2:5 with points palette pt 7 notitle

unset multiplot

