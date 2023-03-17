if [ "$1" = '-h' ];
	then
	echo "order:\n1 - ydiss = 0 gamma; 1 w\n2 - Lsize\n3 - bspacing"
	echo "4 - kappa1\n5 - kappa2\n6 - gamma(y0) or w(y1)"
	exit
fi

if [ ! $6 ]; then
	echo 'ERROR FEW PARAMETERS - 6 MINIMUM'
	echo "write the parameters in this way:"
	echo "order:\n1 - ydiss = [ 0 ---> gamma; 1 ---> w ]\n2 - Lsize\n3 - bspacing"
	echo "4 - kappa1\n5 - kappa2\n6 - gamma(y0) or w(y1)"
	exit
fi

gfortran bkitaev.f90 lapack_tebd.lib.f -o bogol.o
y=$1
L=$2
b=$3
k1=$4
k2=$5
gw=$6

echo "${y}" > init.in #	!ydiss = 0 for gamma & 1 for w
echo "$L" >> init.in #       !Lsize
echo "$b" >> init.in #        !bspacing
echo "$k1" >> init.in #        !kappa1
echo "$k2" >> init.in #        !kappa2
echo "$gw" >> init.in #        !gamma(y0) or w(y1)

./bogol.o #>> run.log
data=$(date)
echo "y$y L$L b$b k1$k1 k2$k2 gw$gw --- END: $data" >> run.log
