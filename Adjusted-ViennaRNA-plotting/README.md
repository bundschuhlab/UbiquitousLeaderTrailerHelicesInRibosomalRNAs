# Adjusted ViennaRNA file for RNAClust mLocarna plots
The "plot_utils.c" file was adjusted for mLocarna plotting to change the greying out of base pairs to be from a numerical number of mismatches to a fraction (full color indicating up to 0.1 fraction mismatches, lightest tint indicating at least 0.3 fraction mismatches). Please replace the plot_utils.c file in "src/ViennaRNA/plotting" prior to building ViennaRNA on your system to replicate our code. 

Note that RNAClust calls mLocarna for final consensus structure building.

ViennaRNA citation:
Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and Hofacker, Ivo L.
ViennaRNA Package 2.0
Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26
