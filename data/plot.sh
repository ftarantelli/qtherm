if [ $2 ]
then
	gnuplot -c plot.plt $1 $2
else
	gnuplot -c plot.plt $1
fi
epstopdf graph.eps
rm graph.eps
