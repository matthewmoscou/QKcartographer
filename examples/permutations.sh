mkdir permutations

for ((  i = 0 ;  i < 1000;  i++  ))
do
	~/bin/Prune -A -i qtlcart.cro -b 2 -t 100 -V
	~/bin/SRmapqtl -A -i qtlcart.crb -t 100 -u 5 -M 2 -V
	~/bin/Zmapqtl -A -i qtlcart.crb -t 1000 -M 6 -V
	mv qtlcart.z permutations/qtlcart_$i.z
done
