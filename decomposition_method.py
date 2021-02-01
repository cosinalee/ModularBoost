import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.cluster
import rpy2.robjects as ro
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from rpy2.robjects.packages import importr


def standardize(X):
    return (X - X.mean()) / (X.std())


def ica_fdr(E, components_num, pic_S=None, pic_A=None, qvalcutoff=None, pic_S_path=None, pic_A_path=None):
    '''
    :param E: Gene expression data, tpye = pandas.DataFrame;
    :param components_num: Number of gene modules, tpye = int;
    :param pic_S: Plot the heat map of the S matrix (p-value), default = False;
    :param pic_A: Plot the heat map of the A matrix (p-value), default = False;
    :param qvalcutoff: Cutoff value, tpye = float;
    :return: Predicted gene module, tpye = list;
    '''

    source = _ica_fastica(E, components_num, pic_S, pic_A, pic_S_path, pic_A_path)
    modules_fdr1 = _ica_fdrtool(E, source, qvalcutoff)

    return modules_fdr1


def ica_fdr_signed(E, components_num, pic_S=None, pic_A=None, qvalcutoff=None, pic_S_path=None, pic_A_path=None):
    source = _ica_fastica(E, components_num, pic_S, pic_A, pic_S_path, pic_A_path)
    modules = _ica_fdrtool_signed(E, source, qvalcutoff)

    return modules


def ica_zscore(E, components_num, pic_S=None, pic_A=None, stdcutoff=None, pic_S_path=None, pic_A_path=None):
    source = _ica_fastica(E, components_num, pic_S, pic_A, pic_S_path, pic_A_path)
    modules = _ica_zscore(E, source, stdcutoff)

    return modules


def pca(E, components_num, qvalcutoff):
    source = _pca(E, components_num)
    modules = _ica_fdrtool(E, source, qvalcutoff)

    return modules


def kmeans(E, k, max_iter=300, n_init=10, seed=None, **kwargs):
    kmeans = sklearn.cluster.KMeans(n_clusters=int(k), max_iter=int(max_iter), n_init=int(n_init), random_state=seed)
    kmeans.fit(standardize(E).T)
    modules = convert_labels2modules(k, kmeans.labels_, E.columns)

    return modules


def _ica_fastica(E, components_num, pic_S=None, pic_A=None, pic_S_path=None, pic_A_path=None):
    '''
    :param E: Gene expression data, tpye = pandas.DataFrame;
    :param components_num: Number of gene modules, tpye = int;
    :param pic_S: Plot the heat map of the S matrix (p-value), default = False;
    :param pic_A: Plot the heat map of the A matrix (p-value), default = False;
    :param pic_S_path: Path of S pic, type = str;
    :param pic_A_path: Path of A pic，type = str;
    :return: S matrix, type = numpy.ndarray;
    '''

    ica = FastICA(n_components=components_num, max_iter=20000, tol=0.001)
    source = ica.fit_transform(standardize(E).T)
    A_ = ica.mixing_

    if pic_S:
        dfsource = pd.DataFrame(source)
        sns.clustermap(dfsource, method='ward', metric='euclidean')
        plt.savefig(pic_S_path, dpi=400)
    if pic_A:
        dfA_ = pd.DataFrame(A_)
        sns.clustermap(dfA_.T, method='ward', metric='euclidean')
        plt.savefig(pic_A_path, dpi=400)

    return source


def _ica_fdrtool(E, source, qvalcutoff):
    '''
    :param E: Gene expression data, tpye = pandas.DataFrame;
    :param source: S matrix, type = numpy.ndarray;
    :param qvalcutoff: Cutoff value, tpye = float;
    :return: Predicted gene module, tpye = list;
    '''
    importr("fdrtool", lib_loc="C:\Program Files\R\R-3.3.3\library")  # 导入R library中的fdr包
    rfdrtool = ro.r["fdrtool"]

    modules = []

    for source_row in source.T:
        rresults = rfdrtool(ro.FloatVector(source_row), plot=False, cutoff_method="fndr", verbose=False)
        qvals = np.array(rresults.rx2("qval"))
        genes = E.columns[qvals < qvalcutoff]
        genes_val = genes.values
        mod_temp = [genes_val[i] for i in range(len(genes_val))]
        modules.append(mod_temp)

    return modules


def _ica_fdrtool_signed(E, source, qvalcutoff):
    importr("fdrtool", lib_loc="C:\Program Files\R\R-3.3.3\library")  # Import the fdr package in the R library
    rfdrtool = ro.r["fdrtool"]

    modules = []

    for source_row in source.T:
        rresults = rfdrtool(ro.FloatVector(source_row), plot=False, cutoff_method="fndr", verbose=False)
        qvals = np.array(rresults.rx2("qval"))

        genes = E.columns[(qvals < qvalcutoff) & (source_row > source_row.mean())]

        genes_val = genes.values
        mod_temp = [genes_val[i] for i in range(len(genes_val))]
        modules.append(mod_temp)

        genes2 = E.columns[(qvals < qvalcutoff) & (source_row < source_row.mean())]

        genes_val2 = genes2.values
        mod_temp2 = [genes_val2[i] for i in range(len(genes_val2))]
        modules.append(mod_temp2)

    return modules


def _ica_zscore(E, source, stdcutoff):
    modules = []
    for source_row in source.T:
        genes = E.columns[source_row < -source_row.std() * stdcutoff].tolist() + E.columns[
            source_row > +source_row.std() * stdcutoff].tolist()
        modules.append(genes)

    return modules


def _pca(E, components_num):
    pca = PCA(n_components=int(components_num))

    source = pca.fit_transform(standardize(E).T)

    return source


def convert_labels2modules(k, labels, G):
    module = []
    names = globals()
    for i in range(k):
        names['temp' + str(i)] = []

    G_value = G.values

    for i in range(len(labels)):
        names['temp' + str(labels[i])].append(G_value[i])
    for i in range(k):
        module.append(names['temp' + str(i)])

    return module
