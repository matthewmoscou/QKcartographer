# QKcartographer
A set of scripts for processing output from QTL Cartographer. The current set of scripts are built around analyzing an F2 population. Adaptation for doubled-haploid and recombinant inbred line populations is current underway.

## Scripts
<i>QKcartographer_preprocess.py</i> converts tabular flat text files into formatted input files for QTL Cartographer.

<i>QKcartographer_permutations.py</i> identifies experiment-wise thresholds based on permuted data sets from QTL Cartographer.

<i>QKcartographer_visualization.py</i> generates figures in PNG or postscript format for publication purposes.

<i>QKcartographer_epistasis.py</i> parses QTL Cartographer output files and permits curation of significant QTLs. Optional command to generate scripts for epistasis analysis with R/qtl.

## Example
Example set of commands to run the set of scripts in Linux:
```bash
python QKcartographer_preprocess.py DxM SF2 DxM_phenotypic_data.txt DxM_genetic_map.txt

mkdir analysis

Rmap -A -i DxM_Rmap.inp -f 2 -V
Rcross -A -i DxM_Rcross.inp -V
SRmapqtl -A -i qtlcart.cro -t 1000 -u 5 -M 2 -V
Zmapqtl -A -i qtlcart.cro -t 1000 -M 6 -V
Eqtl -A -S 8 -H 10 -M 6 -V
mv qtlcart.eqt analysis/qtlcart_H0H1.eqt
Eqtl -A -S 8 -H 30 -M 6 -V
mv qtlcart.eqt analysis/qtlcart_H0H3.eqt
cp qtlcart.cro qtlcart.map qtlcart.sr qtlcart.z analysis/
rm qtlcart.sr qtlcart.z

bash permutations.sh

python QKcartographer_permutations.py DxM 0.95 1000
python QKcartographer_visualization.py DxM SF2 0.95
python QKcartographer_epistasis.py DxM SF2 0.95
```
This example uses a genetic map with cM positions generated using the Kosambi function.

## Supported population type
|Design                            |Code      |Example|Description             |
|:---------------------------------|:--------:|:-----:|:-----------------------|
|Backcross <i>i</i> population     |B<i>i</i> |B1     |Backcross once to mother|
|F<i>i</i> population              |SF<i>i</i>|SF2    |F2 population           |
|Doubled haploid population        |RI0       |RI0    |DH population           |
|Recombinant inbred line population|RI1       |RI1    |RIL population          |

## Permutations
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
