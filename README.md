# QKcartographer
A set of scripts for processing QTL Cartographer output.

QKcartographer_preprocess.py converts tabular flat text files into formatted input files for QTL Cartographer.

QKcartographer_permutations_F2.py identifies experiment-wise thresholds based on permuted data sets from QTL Cartographer.

Permutations for composite interval mapping can be performed using the following bash shell commands:

```bash
for ((  i = 0 ;  i < 1000;  i++  ))
do
  ~/bin/Prune -A -i qtlcart.cro -b 2 -t 100 -V
	~/bin/SRmapqtl -A -i qtlcart.crb -t 100 -u 5 -M 2 -V
	~/bin/Zmapqtl -A -i qtlcart.crb -t 1000 -M 6 -V
	mv qtlcart.z qtlcart_$i.z
done
```

In this case, permutations are performed with reselection, where SRmapqtl is ran on each individual permuted data set (see Lauter et al. (2008) The Plant Genome http://dx.doi.org/10.3835/plantgenome2008.06.0385).

For interval mapping, the bash shell commands are the following:

```bash
for ((  i = 0 ;  i < 1000;  i++  ))
do
  ~/bin/Prune -A -i qtlcart.cro -b 2 -t 100 -V
	~/bin/Zmapqtl -A -i qtlcart.crb -t 1000 -M 3 -V
	mv qtlcart.z qtlcart_$i.z
done
```
