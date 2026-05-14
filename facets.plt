
# ARG0 = this file name
# ARG1 = window number
# ARG2 = input file name
# ARG3 = Title
# ARG4 = number of variables
# ARG5 = 0 => "sigma", 1 => "tau"

if(ARGC < 5) { 
    print("Too few arguments.".ARGC);
    exit;
}

plot_title = ARG3;
n_vars = int(ARG4);
q_idx = n_vars + 2;
w_idx = n_vars + 3;
fn = ARG2;

set palette rgb 33,13,10;
set grid
set tics nomirror

unset tics;
unset border;
unset colorbox;

bw_factor = 50;

set term qt int(ARG1) size 1000,800;

set origin 0,0;
set size 1,1;
set multiplot layout n_vars,n_vars rowsfirst scale 1.0, 1.0 title plot_title;

stats [*:*][*:*] fn using (exp(column(q_idx))) name "Q";
stats [*:*][*:*] fn using (exp(column(w_idx))) name "W";

do for [v=1:n_vars:1] {
    stats [*:*][*:*] fn using (column(v));
    binwidth = (STATS_max - STATS_min) / bw_factor;
    
    do for [r=1:n_vars:1] {
        if(r<v) {

            # Scatter plot
            plot fn using (column(r)):(column(v)):(exp(column(w_idx))) with points palette pt 7 notitle;

        } else if(r==v) {

            # Histograms
            if(r==1)      { set title "a"; }
            else if(r==2) { set title "b"; }
            else if(r==3) { set title "c"; }
            else if(r==4) { set title "d"; }
            else if(r==5) { set title "e"; }
            else if(r==6) { set title "f"; }
            if(r==n_vars) {
                if(ARG5 == 0) { set title "sigma"; }
                else if(ARG5 == 1) { set title "tau"; }
            }

            if(n_vars < 5) {
                plot fn using (bin(column(v), binwidth)):(exp(column(q_idx)) / Q_sum) smooth kdensity bandwidth binwidth with l lw 4 lc rgb "#77FF0000" title "Q", \
                     fn using (bin(column(v), binwidth)):(exp(column(w_idx)) / W_sum) smooth kdensity bandwidth binwidth with l lw 4 lc 19 title "P";
            } else {
                plot fn using (bin(column(v), binwidth)):(exp(column(q_idx)) / Q_sum) smooth kdensity bandwidth binwidth with l lw 4 lc rgb "#77FF0000" notitle, \
                     fn using (bin(column(v), binwidth)):(exp(column(w_idx)) / W_sum) smooth kdensity bandwidth binwidth with l lw 4 lc 19 notitle;
            }

            unset title;

        } else if(r>v) {

            # Correlation
            stats [*:*][*:*] fn using (column(v)):(column(r));
            set label 1 at 0.5,0.5 sprintf("%.4f", STATS_correlation) font "Courier,14";
            plot [0:1] [0:1] NaN title "";
            unset label 1

        }
    }
}

unset multiplot

