# InjectorEmulation
The codes implement the method proposed by Mak, S. et al. (2017+) &lt;arXiv:1611.07911>

To reproduce Figures 6 and 7, and the corresponding values for Figure 8 in the paper, run the file “example2.R”. The R package “InjectorEmulationHelper” needs to be installed from the source package file “InjectorEmulationHelper.tar.gz”. The reproduced results for Figure 6 and 7 are outputted to the directories “output2/POD_expansion/6” and “output2/validation/comparison”, while the results for Figure 8 are printed on the console.  

Because the full simulation data is too large to accompany the manuscript, we have provided a simple example in “example1.R” which uses a simulated toy data to illustrate the full emulation procedure. The raw data can be changed by following the same format in the directory “data1/rawdata folder”: first three columns correspond to the x-, y- and z-coordinates, and other columns are the responses. The experimental design can also be changed from the file “data1/DoE.csv”. The geometry design and raw data for testing validation cases can also be modified in the file “data1/newX.csv” and the folder “data1/testdata”.
