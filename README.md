# Bulk-Sample-Blending
Optimizes a bulk sample of material for which each sample in a set was sub-sampled and assayed thus providing a statistical distribution of the grades. Proportions to take from each sample are optimized to create a blended quantity that will meet desired grade targets.

The main program reads an input file (blend5b.csv) and determins the optimal blend.
The input file is a comma-separated file with the following columns:
Sample source
Sample label
Sample lithology
Sample bin #
# of sample simulatinos
Sample weight (kg)

then for each simulation:
Head iron content (Fe %)
Davis Tube weight recovery (DTWR %)
Davis Tube iron content (FeC %)
Davis Tube silica content (SiC %)

There is one row for each entry consisting of the above format.

