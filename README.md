BioPython scripts elaborated with the help of Sophie Galina &amp; Camille Roux to detect non-synonymous substitutions differently fixed between four different lineages of Silene nutans, in plastid and nuclear genes alignments. 

>**Authors** : Sophie Gallina 1, Camille Roux 1, Zoé Postel 1

>**Affiliations** :
    1. CNRS, Univ. Lille, UMR 8198 – Evo-Eco-Paleo, F-59000 Lille, France

**Corresponding author**: zoe.postel@univ-lille.fr<br />
&nbsp;



# Variant detection for nuclear genes

Script using the output file of Read2snps program (vcf file). 

_Command line to run the script_ : 

./2_select_variants.py --config config.yaml --xls nucl_variants_v2_test.xls

>-- config : give config file for the different direction and file to use

>--xls : give the name of the output file

&nbsp;
&nbsp;


# Variant detection for plastid genes

Directly on aligned and clean multiple sequence alignments in fasta format. 

_Command line to run the script_ : 

./chloroDec_v2.py --samples samples.csv --fasta *.fasta --xls result.xls

>--samples : give a sample file with the lineages they belong to

>--fasta : give the aligned and clean fasta file on which the non-synonymous substitutions will be detected

>--xls : give the name of the result file
