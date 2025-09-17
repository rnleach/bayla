set datafile separator ","

bin(x, width) = width * floor(x / width) + width / 2.0
set style fill solid

# Test Data ----------------------------------------------------------------------------------------------------------------
fn = "test_data.csv"

set key left
set xrange [0.09:3.01]
plot fn using 1:2 with points pt 7 title "Raw Data"
replot log(0.9) + (x - 0.9) / 0.9 - (x - 0.9)**2 / (2 * 0.9**2) + (x -  0.9)**3 / (3 * 0.9**3) lw 4 dt 2 title "True Function"
replot log(x) lw 2 title "Log Model"
replot 7.180454e-01 lw 2 title "Constant"
replot 1.409042e+00 * x + -1.465970e+00 lw 2 title "Linear"
replot 2.848551e-01 * x**2 + 5.259913e-01 * x + -9.914793e-01 lw 2 title "Second Order"
replot 4.245637e-01 * x**3 + -1.689366e+00 * x**2 + 3.023451e+00 * x + -1.700497e+00 lw 2 title "Third Order"

unset xrange

# Simple Model -------------------------------------------------------------------------------------------------------------
fn = "simple_model.csv"
set term qt 1 size 1000,800

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
set term qt 2 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 2,2 rowsfirst scale 0.9, 0.9 title "Log Model Test"

set title "Q"
plot fn using 1:(exp($3)):(exp($3)) with points palette pt 7 notitle
set title "P"
plot fn using 1:(exp($2)):(exp($2)) with points palette pt 7 notitle
set title "Credible Interval"
plot fn using 1:(exp($2)):5 with points palette pt 7 notitle
set title "W"
plot fn using 1:(exp($4)):(exp($4)) with points palette pt 7 notitle

unset multiplot

# Options for future plots -------------------------------------------------------------------------------------------------
unset tics
unset border
unset colorbox

binwidth = 0.02 # Define the width of your bins

# Constant Model -----------------------------------------------------------------------------------------------------------
fn = "constant_model.csv"
set term qt 3 size 1000,800

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
set term qt 4 size 1000,800

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

# Second Order Model -------------------------------------------------------------------------------------------------------
fn = "second_order_model.csv"
set term qt 5 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 4,4 rowsfirst scale 1.0, 1.0 title "Second Order Model Test"

set title "V0"
plot fn using (bin($1, binwidth)):(exp($7)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
set multiplot next
set multiplot next

unset title
plot fn using 1:2:(exp($7)) with points palette pt 7 notitle
set title "b"
plot fn using (bin($2, binwidth)):(exp($7)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
set multiplot next

unset title
plot fn using 1:3:(exp($7)) with points palette pt 7 notitle
plot fn using 2:3:(exp($7)) with points palette pt 7 notitle
set title "c"
plot fn using (bin($3, binwidth)):(exp($7)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next

unset title
plot fn using 1:4:(exp($7)) with points palette pt 7 notitle
plot fn using 2:4:(exp($7)) with points palette pt 7 notitle
plot fn using 3:4:(exp($7)) with points palette pt 7 notitle
set title "sigma"
plot fn using (bin($4, binwidth)):(exp($7)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle

unset multiplot

# Third Order Model --------------------------------------------------------------------------------------------------------
fn = "third_order_model.csv"
set term qt 6 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 5,5 rowsfirst scale 1.0, 1.0 title "Third Order Model Test"

set title "V0"
plot fn using (bin($1, binwidth)):(exp($8)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
set multiplot next
set multiplot next
set multiplot next

unset title
plot fn using 1:2:(exp($8)) with points palette pt 7 notitle
set title "b"
plot fn using (bin($2, binwidth)):(exp($8)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
set multiplot next
set multiplot next

unset title
plot fn using 1:3:(exp($8)) with points palette pt 7 notitle
plot fn using 2:3:(exp($8)) with points palette pt 7 notitle
set title "c"
plot fn using (bin($3, binwidth)):(exp($7)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next
set multiplot next

unset title
plot fn using 1:4:(exp($8)) with points palette pt 7 notitle
plot fn using 2:4:(exp($8)) with points palette pt 7 notitle
plot fn using 3:4:(exp($8)) with points palette pt 7 notitle
set title "d"
plot fn using (bin($4, binwidth)):(exp($8)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle
set multiplot next

unset title
plot fn using 1:5:(exp($8)) with points palette pt 7 notitle
plot fn using 2:5:(exp($8)) with points palette pt 7 notitle
plot fn using 3:5:(exp($8)) with points palette pt 7 notitle
plot fn using 4:5:(exp($8)) with points palette pt 7 notitle
set title "sigma"
plot fn using (bin($5, binwidth)):(exp($8)) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle

unset multiplot

