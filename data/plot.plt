reset

set terminal postscript eps enhanced color dashed font 'Arial, 23'
set output 'graph.eps'


set multiplot


set key samplen 3        # length of line's type
set key font 'Arial, 27'


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
set xlabel "{/Symbol Q}" offset 0, 0.4
set ylabel "L C(x, y)" offset 1.4, 0

set ylabel "L n_{s, eq}" offset 1.4, 0
set xlabel "{/Symbol X}" offset 0, 0.4

namef = "*l"

array length[5] = [11, 30, 60, 90, 120]
array dfile1[5]
array dfile2[5]
array wfile[5]

if(ARG1) {namef = ARG1;}

do for [N=1:5] {
dfile1[N] = sprintf("%s%d.dat", ARG1, length[N]);
}

do for [N=1:5] {
dfile2[N] = sprintf("q3thk-1q0e1.00t1.00S1.000l%d.dat", length[N]);
}

if(ARG2){
do for [N=1:5] {
dfile1[N] = sprintf("%s%d.dat%d", ARG1, length[N], int(ARG2));
}
;}

do for [N=1:2] {
wfile[N+1] = sprintf("L = %d  ", length[N+1]);	#"L = %d  "
wfile[N+3] = sprintf("L = %d  ", length[N+3]);
}
wfile[5] = sprintf("L = %d", length[5]);
wfile[3] = sprintf("L = %d  ", length[3]);

#set key horizontal
set key left top
set ytics 0.2
set xtics 5
#set yrange [-0.4:0.4]
#set xrange [0:0.4]

#set label 1 at graph 0.4, graph 0.96 "{/Arial=30 {/Symbol k}_i = -1}" tc lt 8
set label 2 at graph 0.8, graph 0.06 "{/Arial=30 {/Symbol F} = 0}" tc lt 8
#set label 1 at graph 0.8, graph 0.16 "{/Arial=30 {/Symbol F}_i = -1}" tc lt 8
#set label 3 at graph 0.3, graph 0.06 "{/Arial=30 {/Symbol g} = 1}" tc lt 8
#set label 3 at graph 0.65, graph 0.95 "{/Arial=30 {/Symbol G} = 10}" tc lt 3
#set label 5 at graph 0.45, graph 0.78 "{/Arial=30 {/Symbol G} = 1}" tc lt 8
#set label 4 at graph 0.8, graph 0.26 "{/Arial=30 {/Symbol X} = 1}" tc lt 8
#set label 5 at graph 0.8, graph 0.86 "{/Arial=30 {/Symbol K_b = 2}" tc lt 8

#set arrow from graph 0, first 0.026873444 to graph 1, first 0.026873444 nohead lt rgb 'dark-green' lw 4 dt '_'

set errorbars

plot    dfile1[5] using ($1*length[5]**(0.)):($4*length[5]**(0.)) t wfile[5] w l lw 4. lt 8 dt 1,\
        dfile1[4] using ($1*length[4]**(0.)):($4*length[4]**(0.)) t wfile[4] w l lw 4. lt 1 dt 7,\
        dfile1[3] using ($1*length[3]**(0.)):($4*length[3]**(0.)) t wfile[3] w l lw 4. lt rgb 'red' dt 6,\
        dfile1[2] using ($1*length[2]**(0.)):($4*length[2]**(0.)) t wfile[2] w l lw 4. lt rgb 'blue' dt 4,\
        dfile1[1] using ($1*length[1]**(0.)):($4*length[1]**(0.)) t wfile[1] w l lw 4. lt rgb 'dark-green' dt 5,\
        #dfile2[5] using ($1*length[5]**(0.)):($4*length[5]**(0.)) t 'L = 120' w l lw 4. lt 3 dt 1,\
        #dfile2[4] using ($1*length[4]**(0.)):($4*length[4]**(0.)) t '' w l lw 4. lt 1 dt 7,\
        #dfile2[3] using ($1*length[3]**(0.)):($4*length[3]**(0.)) t '' w l lw 4. lt rgb 'red' dt 6,\
        #dfile2[2] using ($1*length[2]**(0.)):($4*length[2]**(0.)) t '' w l lw 4. lt rgb 'blue' dt 4,\
        #dfile2[1] using ($1*length[1]**(0.)):($4*length[1]**(0.)) t '' w l lw 4. lt rgb 'dark-green' dt 5,\

        #'qthk-1q1e0l60.dat' using ($1/60):($4*60) t 'T = 0' w l lw 4. lt 1 dt 2,\
        #dfile1[4] using ($4):($4*length[4]**(0.)) t wfile[4] w l lw 4. lt 8 dt 1,\

unset margin
set size 0.525,0.4
set origin 0.45,0.58

unset xlabel
unset ylabel
unset title

set ytics font "Arial,23"
set xtics  font "Arial,23"
#set xtics 0.02#, _, _  # offset _, _
unset key
#set xrange
unset label
unset xrange
unset yrange

set ytics 0.1
set xtics 0.01

set xlabel font "Arial,26"
#set ylabel font "Arial,27"
set label font "Arial,25"

set ytics font "Arial,20"
set xtics  font "Arial,20"

#set format x "%.3f"
#set yrange [0.76:0.]
set xrange [0.00:0.018]
set xlabel "1/L" offset 6, 1.8
#set ylabel "D({/Arial=9 L/2, t}) - D({/Arial=9 L/2, 0})" offset 3.5,0
#set xlabel "t L^{- 3}" offset 0,0.8
#set ylabel "L < C^+_{3L/8} C_{5L/8} + h.c. >" offset 2.5,0

set label 1 at screen 0.62, screen 0.89 "{/Symbol=28 Q} = 0.5"

inset = sprintf("inset%s.dat", ARG1);

#f(x)=a*x + b
#fit f(x) inset u (1./$1):($5*$1**(0.)) via a,b


#plot  f(x) notitle w l dashtype '_' lt 1, inset using (1./$1):($5*$1**(0.)) t '' pt 7 ps 2 # lt 8

unset multiplot
