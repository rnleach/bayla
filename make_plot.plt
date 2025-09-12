set datafile separator ","

bin(x, width) = width * floor(x / width) + width / 2.0
set style fill solid

# Simple Model -------------------------------------------------------------------------------------------------------------
fn = "simple_model.csv"
set term qt 0 size 1000,800

set palette rgb 33,13,10
set grid
set tics nomirror

set origin 0,0
set size 1,1
set multiplot layout 2,2 rowsfirst scale 0.9, 0.9 title "Simple Model Test"

set title "Q"
plot fn using 1:2:(exp($4)) with points palette pt 7 notitle
set title "P"
plot fn using 1:2:(exp($3)) with points palette pt 7 notitle
set title "Credible Interval"
plot fn using 1:2:6 with points palette pt 7 notitle
set title "W"
plot fn using 1:2:(exp($5)) with points palette pt 7 notitle

unset multiplot

# Log Model ----------------------------------------------------------------------------------------------------------------
fn = "log_model.csv"
set term qt 1 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 2,2 rowsfirst scale 0.9, 0.9 title "Log Model Test"

set title "Q"
plot fn using 1:3:(exp($3)) with points palette pt 7 notitle
set title "P"
plot fn using 1:2:(exp($2)) with points palette pt 7 notitle
set title "Credible Interval"
plot fn using 1:2:5 with points palette pt 7 notitle
set title "W"
plot fn using 1:4:(exp($4)) with points palette pt 7 notitle

unset multiplot

# Options for future plots -------------------------------------------------------------------------------------------------
unset tics
unset border
unset colorbox

binwidth = 0.02 # Define the width of your bins

# Constant Model -----------------------------------------------------------------------------------------------------------
fn = "constant_model.csv"
set term qt 2 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 2,2 rowsfirst scale 0.9, 0.9 title "Constant Model Test"

set title "V0"
plot fn using (bin($1, binwidth)):(exp($5)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
unset title

plot fn using 1:2:(exp($5)) with points palette pt 7 notitle
set title "sigma"
plot fn using (bin($2, binwidth)):(exp($5)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle

unset multiplot

# Linear Model -------------------------------------------------------------------------------------------------------------
fn = "linear_model.csv"
set term qt 3 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 3,3 rowsfirst scale 1.0, 1.0 title "Linear Model Test"

set title "V0"
plot fn using (bin($1, binwidth)):(exp($6)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
set multiplot next

unset title
plot fn using 1:2:(exp($6)) with points palette pt 7 notitle
set title "b"
plot fn using (bin($2, binwidth)):(exp($6)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next

unset title
plot fn using 1:3:(exp($6)) with points palette pt 7 notitle
plot fn using 2:3:(exp($6)) with points palette pt 7 notitle
set title "sigma"
plot fn using (bin($3, binwidth)):(exp($6)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle

unset multiplot

