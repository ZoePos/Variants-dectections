============== Variant detection for nuclear genes =======================

Script using the output file of Read2snps program (vcf file). 

Command line to run the script : 
./2_select_variants.py --config config.yaml --xls=nucl_variants_v2_test.xls

-- config : give config file for the different direction and file to use
--xls : give the name of the output file


============== Variant detection for plastid genes =======================

Directly on aligned and clean multiple sequence alignments in fasta format. 

Command line to run the script : 
./chloroDec_v2.py --samples samples.csv --fasta *.fasta --xls result.xls

--samples : give a sample file with the lineages they belong to
--fasta : give the aligned and clean fasta file on which the non-synonymous substitutions will be detected
--xls : give the name of the result file
