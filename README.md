# QKcartographer
A set of scripts for processing QTL Cartographer output.

<i>QKcartographer_preprocess.py</i> converts tabular flat text files into formatted input files for QTL Cartographer.

<i>QKcartographer_permutations_F2.py</i> identifies experiment-wise thresholds based on permuted data sets from QTL Cartographer.

<i>QKcartographer_visualization.py</i> generates figures in PNG or postscript format for publication purposes.

<i>QKcartographer_epistasis.py</i> parses QTL Cartographer output files and permits curation of significant QTLs. Optional command to generate scripts for epistasis analysis with R/qtl.

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
python QKcartographer_visualization.py SF2 0.95
```
