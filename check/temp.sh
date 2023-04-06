for ii in `seq 1 1 4`
do
	size=$(echo "${ii}*50" |bc -l)
	Temper="0.5"
	file="data/equil0T${Temper}L${size}.dat"
	#./therm.out k0. T${Temper} L${size} y5 > ${file} &
	eta="1"
	file="data/equil0e${eta}L${size}.dat"
	./therm.out k0. e${eta} L${size} y5 > ${file} &
done
