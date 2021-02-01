# ModularBoost
ModularBoost: An efficient network inference algorithm based on Module Decomposition
## Explanation for the datasets
Cruted, PIDC, SCODE, and DREAM5 datasets come from the previous work of researchers, which can be trace back to Reference. It should be noted that the 'Goldmodule.txt' in the Curated and SCODE datasets are gold standard of module detection. They were extracted by graph theory or community detection methods.
## Explanation for the code
We provide code for gene module detection, including ICA-FDR, ICA-FDR2, ICA-zscore, PCA and K-means.  
The code developed and tested in [Python 3.6](https://www.python.org/downloads/release/python-360/) and requires:  
[numpy](https://numpy.org/), for data structures;  
[pandas](https://pandas.pydata.org/), for managing data;  
[seaborn] (http://seaborn.pydata.org/), for heatmaps of S and A matrix;  
[sklearn](https://scikit-learn.org/stable/), for decomposition algorithm;  
[rpy2](https://pypi.org/project/rpy2/), for FDR functions;  
[matplotlib](https://matplotlib.org/), for plotting;
