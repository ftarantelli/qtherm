reset

set terminal postscript eps enhanced color dashed font 'Arial, 23'
set output 'graph.eps'


set multiplot


set key samplen 3        # length of line's type
set key font 'Arial, 27'

#set logscale x 10
set logscale y 10

#set format x "%2.0t{/Symbol \264}10^{%L}"
#set format y "%2.4t{/Symbol \264}10^{%L}"
set format y "10^{%L}"

set lmargin at screen 0.14
set bmargin at screen 0.129
set rmargin at screen 0.965
set mxtics 4
set xtics scale 2.3, 1.
set mytics 4
set ytics scale 2.3, 1.

set xtics add ("0" 0)
set ytics add ("0" 0)

set ytics font "Arial,20"
set xtics  font "Arial,20"

set ylabel font 'Arial, 35'
set xlabel font 'Arial, 35'

set label font "Arial,30"

set xlabel "{/Symbol m - m}_c L^{y_m}" offset 0, 0.8
set xlabel "t" offset 0, 0.3
set ylabel "C(x, y)" offset 2.4, -0.2
set ylabel "n_s - n_{s, eq}(T)" offset 1.8, 0

namef = "*l"

array length[5] = [0., 3, 2, 1, 0]
array dfile1[5]
array wfile[5]

if(ARG1) {namef = ARG1;}

do for [N=1:5] {
dfile1[N] = sprintf("%s%d.dat", ARG1, length[N]);
}

if(ARG2){
do for [N=1:5] {
dfile1[N] = sprintf("%s%d.dat%d", ARG1, length[N], int(ARG2));
}
;}

do for [N=1:2] {
wfile[N+1] = sprintf("{/Symbol g} = 10^{-%d}", length[N+1]);	#"L = %d  "
wfile[N+3] = sprintf("{/Symbol g} = 10^{-%d}", length[N+3]);
}
wfile[1] = sprintf("{/Symbol g} = 0    ");
wfile[5] = sprintf("{/Symbol g} = 1    ");

#set key horizontal
set key at graph 0.99, graph 0.65
set ytics 0.00001
set xtics 50
set xrange [0:100]
#set yrange [0:0.14]

#set label 1 at graph 0.4, graph 0.96 "{/Arial=30 {/Symbol k}_i = -1}" tc lt 8
set label 2 at graph 0.25, graph 0.16 "{/Arial=30 w = 0}" tc lt 8
set label 1 at graph 0.25, graph 0.06 "{/Arial=30 w_i = -0.01}" tc lt 8
set label 3 at graph 0.25, graph 0.26 "{/Arial=30 L = 60}" tc lt 8
#set label 3 at graph 0.8, graph 0.66 "{/Arial=30 {/Symbol G} = 1 }" tc lt 8
set label 4 at graph 0.8, graph 0.74 "{/Arial=30 T = 2}" tc lt 8
#set label 5 at graph 0.8, graph 0.66 "{/Arial=30 {/Symbol h}_b = 1}" tc lt 8

set errorbars

plot    "q3thk-1q0e0.00t1.00g1.000l60.dat" using ($1*60.**(1.)):(($4-0.027607198)*60**(-1.)) t wfile[5] w l lw 6. lt 8 dt 1,\
        "q3thk-1q0e0.00t1.00g0.100l60.dat" using ($1*60.**(1.)*1):(($4-0.027607198)*60**(-1.)) t wfile[4] w l lw 4. lt 1 dt 7,\
        "q3thk-1q0e0.00t1.00g0.010l60.dat" using ($1*60.**(1.)*1):(($4-0.027607198)*60**(-1.)) t wfile[3] w l lw 4. lt rgb 'red' dt 6,\
        "q3thk-1q0e0.00t1.00g0.001l60.ddat" using ($1*60.**(1.)*1):(($4-0.027607198)*60**(-1.)) t wfile[2] w l lw 4. lt rgb 'blue' dt 4,\
        "q3thk-1q0e0.00t1.00g0.000l60.dat" using ($1*60.**(1.)):(($4-0.027607198)*60**(-1.)) t wfile[1] w l lw 2. lt rgb 'dark-green' dt 2,\

# 0.00046011996666666666 = n_{s, asymptotic}

        #'qthk-1q1e0l60.dat' using ($1/60):(($4-0.027607198)*60) t 'T = 0' w l lw 4. lt 1 dt 2,\
        #dfile1[4] using ($2):($2*length[4]**(0.)) t wfile[4] w l lw 4. lt 8 dt 1,\
unset multiplot
