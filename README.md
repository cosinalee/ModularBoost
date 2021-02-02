# ModularBoost
ModularBoost: An efficient network inference algorithm based on Module Decomposition
## Explanation for the datasets
Cruted[1], PIDC[2], SCODE[3], and DREAM5[4] datasets come from the previous work of researchers, which can be trace back to Reference. It should be noted that the 'Goldmodule.txt' in the Curated and SCODE datasets are gold standard of module detection. They were extracted by graph theory or community detection methods.
## Explanation for the codes
We provide code for gene module detection, including ICA-FDR, ICA-FDR2, ICA-zscore, PCA and K-means.  
The code developed and tested in [Python 3.6](https://www.python.org/downloads/release/python-360/) and requires:  
[numpy](https://numpy.org/), for data structures;  
[pandas](https://pandas.pydata.org/), for managing data;  
[seaborn](http://seaborn.pydata.org/), for heatmaps of S and A matrix;  
[sklearn](https://scikit-learn.org/stable/), for decomposition algorithm;  
[rpy2](https://pypi.org/project/rpy2/), for FDR functions;  
[matplotlib](https://matplotlib.org/), for plotting;
### Reference
[1] Pratapa, Aditya, et al. [Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data.](https://www.nature.com/articles/s41592-019-0690-6)[J] Nature methods 17.2 (2020): 147-154.  
[2] Chan T E, Stumpf M P H, Babtie A C. [Gene regulatory network inference from single-cell data using multivariate information measures](https://www.sciencedirect.com/science/article/pii/S2405471217303861)[J]. Cell systems, 2017, 5(3): 251-267. e3.  
[3] Matsumoto H, Kiryu H, Furusawa C, et al. [SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation](https://academic.oup.com/bioinformatics/article/33/15/2314/3100331?login=true)[J]. Bioinformatics, 2017, 33(15): 2314-2321.  
[4] Marbach D, Costello J C, Küffner R, et al. [Wisdom of crowds for robust gene network inference](https://www.nature.com/articles/nmeth.2016)[J]. Nature methods, 2012, 9(8): 796-804.  
[5] Saelens W, Cannoodt R, Saeys Y. [A comprehensive evaluation of module detection methods for gene expression data](https://www.nature.com/articles/s41467-018-03424-4)[J]. Nature communications, 2018, 9(1): 1-12.
