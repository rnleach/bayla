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
call "facets.plt" 3 "constant_model.csv"                  "Constant Model Test (LaPlace Approx)"           2 0;
call "facets.plt" 4 "constant_tau_model.csv"              "Constant Tau Model Test (LaPlace Approx)"       2 1;
call "facets.plt" 5 "stud_t_constant_model.csv"           "Constant Model Test (Student's-T Approx)"       2 0;
call "facets.plt" 6 "stud_t_constant_tau_model.csv"       "Constant Tau Model Test (Student's-T Approx)"   2 1;

# Linear Model -------------------------------------------------------------------------------------------------------------
call "facets.plt" 7 "linear_model.csv"                    "Linear Model Test (LaPlace Approx)"             3 0;
call "facets.plt" 8 "stud_t_linear_model.csv"             "Linear Model Test (Student's-T Approx)"         3 0;
call "facets.plt" 9 "linear_tau_model.csv"                "Linear Tau Model Test (LaPlace Approx)"         3 1;
call "facets.plt" 10 "stud_t_linear_tau_model.csv"        "Linear Tau Model Test (Student's-T Approx)"     3 1;

# Second Order Model -------------------------------------------------------------------------------------------------------
call "facets.plt" 11 "second_order_model.csv"             "Second Order Model Test (LaPlace Approx)"       4 0;
call "facets.plt" 12 "stud_t_second_order_model.csv"      "Second Order Model Test (Student-T Approx)"     4 0;
call "facets.plt" 13 "second_order_tau_model.csv"         "Second Order Tau Model Test (LaPlace Approx)"   4 1;
call "facets.plt" 14 "stud_t_second_order_tau_model.csv"  "Second Order Tau Model Test (Student-T Approx)" 4 1;

# Third Order Model --------------------------------------------------------------------------------------------------------
call "facets.plt" 15 "third_order_model.csv"             "Third Order Model Test (LaPlace Approx)"         5 0;
call "facets.plt" 16 "stud_t_third_order_model.csv"      "Third Order Model Test (Student-T Approx)"       5 0;
call "facets.plt" 17 "third_order_tau_model.csv"         "Third Order Tau Model Test (LaPlace Approx)"     5 1;
call "facets.plt" 18 "stud_t_third_order_tau_model.csv"  "Third Order Tau Model Test (Student-T Approx)"   5 1;

# Fourth Order Model --------------------------------------------------------------------------------------------------------
call "facets.plt" 19 "fourth_order_model.csv"            "Fourth Order Model Test (LaPlace Approx)"        6 0;
call "facets.plt" 20 "stud_t_fourth_order_model.csv"     "Fourth Order Model Test (Student-T Approx)"      6 0;
call "facets.plt" 21 "fourth_order_tau_model.csv"        "Fourth Order Tau Model Test (LaPlace Approx)"    6 1;
call "facets.plt" 22 "stud_t_fourth_order_tau_model.csv" "Fourth Order Tau Model Test (Student-T Approx)"  6 1;

# Fifth Order Model --------------------------------------------------------------------------------------------------------
call "facets.plt" 23 "fifth_order_model.csv"             "Fifth Order Model Test (LaPlace Approx)"         7 0;
call "facets.plt" 24 "stud_t_fifth_order_model.csv"      "Fifth Order Model Test (Student-T Approx)"       7 0;
call "facets.plt" 25 "fifth_order_tau_model.csv"         "Fifth Order Tau Model Test (LaPlace Approx)"     7 1;
call "facets.plt" 26 "stud_t_fifth_order_tau_model.csv"  "Fifth Order Tau Model Test (Student-T Approx)"   7 1;

