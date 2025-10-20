set datafile separator ","

bin(x, width) = width * floor(x / width) + width / 2.0
set style fill solid

# Test Data ----------------------------------------------------------------------------------------------------------------
fn = "test_data.csv"
set term qt 0 size 1600,800

load "obs_models.plt"

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

# Constant Model -----------------------------------------------------------------------------------------------------------
fn = "constant_model.csv"
set term qt 3 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 2,2 rowsfirst scale 0.9, 0.9 title "Constant Model Test"

set title "V0"
fn2 = "constant_model_marginal_v0.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 1:2
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""
unset label 1

plot fn using 1:2:(exp($5)) with points palette pt 7 notitle
set title "sigma"
fn2 = "constant_model_marginal_sigma.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle

unset multiplot

# Linear Model -------------------------------------------------------------------------------------------------------------
fn = "linear_model.csv"
set term qt 4 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 3,3 rowsfirst scale 1.0, 1.0 title "Linear Model Test"

set title "V0"
fn2 = "linear_model_marginal_v0.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 1:2
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""
unset label 1

plot fn using 1:2:(exp($6)) with points palette pt 7 notitle
set title "b"
fn2 = "linear_model_marginal_b.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 2:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""
unset label 1

plot fn using 1:3:(exp($6)) with points palette pt 7 notitle
plot fn using 2:3:(exp($6)) with points palette pt 7 notitle
set title "sigma"
fn2 = "linear_model_marginal_sigma.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle

unset multiplot

# Second Order Model -------------------------------------------------------------------------------------------------------
fn = "second_order_model.csv"
set term qt 5 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 4,4 rowsfirst scale 1.0, 1.0 title "Second Order Model Test"

set title "V0"
fn2 = "second_order_model_marginal_v0.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 1:2
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

unset title
plot fn using 1:2:(exp($7)) with points palette pt 7 notitle
set title "b"
fn2 = "second_order_model_marginal_b.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 2:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 2:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:3:(exp($7)) with points palette pt 7 notitle
plot fn using 2:3:(exp($7)) with points palette pt 7 notitle
set title "c"
fn2 = "second_order_model_marginal_c.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 3:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

unset title
plot fn using 1:4:(exp($7)) with points palette pt 7 notitle
plot fn using 2:4:(exp($7)) with points palette pt 7 notitle
plot fn using 3:4:(exp($7)) with points palette pt 7 notitle
set title "sigma"
fn2 = "second_order_model_marginal_sigma.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle

unset multiplot

# Third Order Model --------------------------------------------------------------------------------------------------------
fn = "third_order_model.csv"
set term qt 6 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 5,5 rowsfirst scale 1.0, 1.0 title "Third Order Model Test"

set title "V0"
fn2 = "third_order_model_marginal_v0.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 1:2
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:2:(exp($8)) with points palette pt 7 notitle
set title "b"
fn2 = "third_order_model_marginal_b.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 2:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 2:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 2:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:3:(exp($8)) with points palette pt 7 notitle
plot fn using 2:3:(exp($8)) with points palette pt 7 notitle
set title "c"
fn2 = "third_order_model_marginal_c.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 3:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 3:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:4:(exp($8)) with points palette pt 7 notitle
plot fn using 2:4:(exp($8)) with points palette pt 7 notitle
plot fn using 3:4:(exp($8)) with points palette pt 7 notitle
set title "d"
fn2 = "third_order_model_marginal_d.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 4:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:5:(exp($8)) with points palette pt 7 notitle
plot fn using 2:5:(exp($8)) with points palette pt 7 notitle
plot fn using 3:5:(exp($8)) with points palette pt 7 notitle
plot fn using 4:5:(exp($8)) with points palette pt 7 notitle
set title "sigma"
fn2 = "third_order_model_marginal_sigma.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle

unset multiplot

# Fourth Order Model --------------------------------------------------------------------------------------------------------
fn = "fourth_order_model.csv"
set term qt 7 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 6,6 rowsfirst scale 1.0, 1.0 title "Fourth Order Model Test"

set title "V0"
fn2 = "fourth_order_model_marginal_v0.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 1:2
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 1:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:2:(exp($9)) with points palette pt 7 notitle
set title "b"
fn2 = "fourth_order_model_marginal_b.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 2:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 2:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 2:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 2:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:3:(exp($9)) with points palette pt 7 notitle
plot fn using 2:3:(exp($9)) with points palette pt 7 notitle
set title "c"
fn2 = "fourth_order_model_marginal_c.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 3:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 3:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 3:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:4:(exp($9)) with points palette pt 7 notitle
plot fn using 2:4:(exp($9)) with points palette pt 7 notitle
plot fn using 3:4:(exp($9)) with points palette pt 7 notitle
set title "d"
fn2 = "fourth_order_model_marginal_d.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 4:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

stats fn using 4:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:5:(exp($9)) with points palette pt 7 notitle
plot fn using 2:5:(exp($9)) with points palette pt 7 notitle
plot fn using 3:5:(exp($9)) with points palette pt 7 notitle
plot fn using 4:5:(exp($9)) with points palette pt 7 notitle
set title "e"
fn2 = "fourth_order_model_marginal_e.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 5:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:6:(exp($9)) with points palette pt 7 notitle
plot fn using 2:6:(exp($9)) with points palette pt 7 notitle
plot fn using 3:6:(exp($9)) with points palette pt 7 notitle
plot fn using 4:6:(exp($9)) with points palette pt 7 notitle
plot fn using 5:6:(exp($9)) with points palette pt 7 notitle
set title "sigma"
fn2 = "fourth_order_model_marginal_sigma.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle

unset multiplot

# Fifth Order Model --------------------------------------------------------------------------------------------------------
fn = "fifth_order_model.csv"
set term qt 8 size 1000,800

set origin 0,0
set size 1,1
set multiplot layout 7,7 rowsfirst scale 1.0, 1.0 title "Fifth Order Model Test"

set title "V0"
fn2 = "fifth_order_model_marginal_v0.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 1:2
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 1:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 1:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 1:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 1:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 1:7
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:2:(exp($10)) with points palette pt 7 notitle
set title "b"
fn2 = "fifth_order_model_marginal_b.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 2:3
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 2:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 2:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 2:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 2:7
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:3:(exp($10)) with points palette pt 7 notitle
plot fn using 2:3:(exp($10)) with points palette pt 7 notitle
set title "c"
fn2 = "fifth_order_model_marginal_c.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 3:4
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 3:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 3:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 3:7
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:4:(exp($10)) with points palette pt 7 notitle
plot fn using 2:4:(exp($10)) with points palette pt 7 notitle
plot fn using 3:4:(exp($10)) with points palette pt 7 notitle
set title "d"
fn2 = "fifth_order_model_marginal_d.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 4:5
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 4:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 4:7
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:5:(exp($10)) with points palette pt 7 notitle
plot fn using 2:5:(exp($10)) with points palette pt 7 notitle
plot fn using 3:5:(exp($10)) with points palette pt 7 notitle
plot fn using 4:5:(exp($10)) with points palette pt 7 notitle
set title "e"
fn2 = "fifth_order_model_marginal_e.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 5:6
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

stats fn using 5:7
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:6:(exp($10)) with points palette pt 7 notitle
plot fn using 2:6:(exp($10)) with points palette pt 7 notitle
plot fn using 3:6:(exp($10)) with points palette pt 7 notitle
plot fn using 4:6:(exp($10)) with points palette pt 7 notitle
plot fn using 5:6:(exp($10)) with points palette pt 7 notitle
set title "f"
fn2 = "fifth_order_model_marginal_f.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle
unset title

stats fn using 6:7
set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,12"
plot [0:1] [0:1] NaN title ""

unset label 1

plot fn using 1:7:(exp($10)) with points palette pt 7 notitle
plot fn using 2:7:(exp($10)) with points palette pt 7 notitle
plot fn using 3:7:(exp($10)) with points palette pt 7 notitle
plot fn using 4:7:(exp($10)) with points palette pt 7 notitle
plot fn using 5:7:(exp($10)) with points palette pt 7 notitle
plot fn using 6:7:(exp($10)) with points palette pt 7 notitle
set title "sigma"
fn2 = "fifth_order_model_marginal_sigma.csv"
plot fn2 using 1:2 with boxes lw 2 lc 19 notitle

unset multiplot

