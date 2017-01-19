# QKcartographer
A set of scripts for processing QTL Cartographer output.

QKcartographer_preprocess.py converts tabular flat text files into formatted input files for QTL Cartographer.

QKcartographer_permutations_F2.py identifies experiment-wise thresholds based on permuted data sets from QTL Cartographer.

QKcartographer_visualization.py generates figures in PNG or postscript format for publication purposes.

QKcartographer_epistasis.py parses QTL Cartographer output files and permits curation of significant QTLs. Optional command to generate scripts for epistasis analysis with Rqtl.

Permutations for composite interval mapping can be performed using the following bash shell commands:

```bash
mkdir permutations

for ((  i = 0 ;  i < 1000;  i++  ))
do
  ~/bin/Prune -A -i qtlcart.cro -b 2 -t 100 -V
	~/bin/SRmapqtl -A -i qtlcart.crb -t 100 -u 5 -M 2 -V
	~/bin/Zmapqtl -A -i qtlcart.crb -t 1000 -M 6 -V
	mv qtlcart.z permutations/qtlcart_$i.z
done
```

In this case, permutations are performed with reselection, where SRmapqtl is ran on each individual permuted data set (see Lauter et al. (2008) The Plant Genome http://dx.doi.org/10.3835/plantgenome2008.06.0385).

For interval mapping, the bash shell commands are the following:

```bash
mkdir permutations

for ((  i = 0 ;  i < 1000;  i++  ))
do
  ~/bin/Prune -A -i qtlcart.cro -b 2 -t 100 -V
	~/bin/Zmapqtl -A -i qtlcart.crb -t 1000 -M 3 -V
	mv qtlcart.z permutations/qtlcart_$i.z
done
```

Example set of commands to run the set of scripts:
```bash
python QKcartographer_preprocess.py DxM SF2 DxM_phenotypic_data.txt DxM_genetic_map.txt
bash permutation.sh
python QKcartographer_permutations_F2.py 0.95 1000
python QKcartographer_visualization.py DxM SF2 0.95
```
