# Adjusted ViennaRNA file for RNAClust mLocarna plots
The "plot_utils.c" file was adjusted for mLocarna plotting to change the greying out of base pairs to be from a numerical number of mismatches to a fraction (full color indicating up to 0.1 fraction mismatches, lightest tint indicating at least 0.3 fraction mismatches). Please replace the plot_utils.c file in "src/ViennaRNA/plotting" prior to building ViennaRNA on your system to replicate our code. 

Note that RNAClust calls mLocarna for final consensus structure building.