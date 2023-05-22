for script in q3thk-1q0e1.00t1.00S10.000l*
do
	num=$(wc -l < $script)
	num=$(echo "${num}*1/20+4" | bc -l)
	num=${num%.*}
	size=${script##*l}
	size=${size%%.dat}
	if [ "$1" != '-n' ];
	then
		echo -e "\033[35;1m$script\033[0m"
	fi
	values=$(sed -n "${num}p" $script)
	echo -e "$size\t$values"
	#echo -e "\n"
done
