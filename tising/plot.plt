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

set xlabel "h L^{15/8}" offset 0, 0.4
set ylabel "L^{1/8} M_x" offset 3.2, 0


namef = "*l"

array length[4] = [5, 8, 9, 10]
array dfile1[4]
array wfile[4]

if(ARG1) {namef = ARG1;}

do for [N=1:4] {
dfile1[N] = sprintf("%d.dat", length[N]);
}

if(ARG2){
do for [N=1:4] {
dfile1[N] = sprintf("%s%d.dat%d", ARG1, length[N], int(ARG2));
}
;}

do for [N=1:2] {
wfile[N] = sprintf("L = %d", length[N]);	#"L = %d  "
wfile[N+2] = sprintf("L = %d", length[N+2]);
}

set key top left
set ytics 0.03
set xtics 0.1
#set xrange [-5:15]

set label 1 at graph 0.4, graph 0.96 "{/Arial=30 g = 1}" tc lt 8
set label 2 at graph 0.6, graph 0.96 "{/Arial=30 TL = 100}" tc lt 8
#set label 3 at graph 0.4, graph 0.86 "{/Arial=30 hL^{y_h} = 0.01}" tc lt 8

set errorbars

plot dfile1[4] using ($3/length[4]**(0.)):($5*length[4]**(0.)) t wfile[4] w l lw 4. lt 8 dt 1,\
        dfile1[3] using ($3/length[3]**(0.)):($5*length[3]**(0.)) t wfile[3] w l lw 4. lt rgb 'red' dt 6,\
        dfile1[2] using ($3/length[2]**(0.)):($5*length[2]**(0.)) t wfile[2] w l lw 4. lt rgb 'blue' dt 4,\
        dfile1[1] using ($3/length[1]**(0.)):($5*length[1]**(0.)) t wfile[1] w l lw 4. lt rgb 'dark-green' dt 5,\

	 #'qthk-1q1e0l60.dat' using ($3/60):($5*60) t 'T = 0' w l lw 4. lt 1 dt 2,\
     #dfile1[4] using ($3):($4*length[4]**(0.)) t wfile[4] w l lw 4. lt 8 dt 1,\

unset multiplot
