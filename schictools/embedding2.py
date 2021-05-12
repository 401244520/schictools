# -*- coding: utf-8 -*-
# 本模块提供染色体的交互矩阵的图嵌入表示

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
import warnings,time

def fast_oe(x):
    l = len(x)
    y = np.zeros((l,l))
    for i in range(l):
        s = 0
        for j in range(l-i):
            s += (x[j,j+i])
        exc = s/(l-i)
        if exc > 0 :
            for k in range(l-i):
                if x[k,k+i] > 0 :
                    y[k,k+i] = x[k,k+i]/exc
                    y[k+i,k] = x[k,k+i]/exc
    return(y)

def decomposition_mainvector(mat,pc = 1):
    warnings.filterwarnings("ignore") 
    corr = np.corrcoef(fast_oe(mat))  
    corr[np.isnan(corr)] = 0
    pca = PCA(2).fit_transform(corr)
    #plt.fill_between(range(len(pca)),0,pca[:,pc])
    return(pca[:,pc])


def call_compartments2(compartmentpath= "./data/ES-E14_bulk.cis.vecs.tsv",resolution = 1000000):
    # TODO：use cooltools api rather than cli to calculate comparment engivectors
    data = pd.read_table(compartmentpath)
    forward = 1
    new = []
    compartmentname = 0
    for index,values in data.iterrows():
        if forward * values.E2 > 0 :
            new.append(compartmentname)
            forward = values.E2
        elif forward * values.E2 <= 0 :
            compartmentname += 1
            new.append(compartmentname)
            forward = values.E2
        else :
            forward = 0
            new.append(np.nan)
        print(compartmentname)
        len(new)
    data["counts"] = new
    data["start"] = (data["start"]/resolution).astype(int)
    data = data.dropna()
    data[["chrom","start","counts"]].to_csv("{}_2.csv".format(compartmentpath),index=0)
    return(data)

def load_bins(binspath = "./data/data/bins_1mb.csv",resolution = 1000000):
    bins = pd.read_csv(binspath)
    bins["start"] = (bins["start"]/resolution).astype(int)
    bins["counts"] = bins.index
    bins_index = bins.set_index(["chrom","start"])["counts"].to_dict()
    bins.to_csv("./data/bins_1mb.csv",index = 0)


#计算交互矩阵,以bin或compartment为单位
def compartment_concentrate(contacts,comp,resolution = 1000000,weighted = False):
    contacts = contacts.copy()
    compartmentnames = comp.set_index(["chrom","start"])["counts"].to_dict()
    contacts[["pos1","pos2"]] = (contacts[["pos1","pos2"]]/resolution).astype("int")
    matlen= int(max(compartmentnames.values()))
    mat = np.zeros((matlen,matlen))
    #TODO:this step is time consuming, need to find better algorism
    for index,row in contacts.iterrows():
        pos1,pos2 = compartmentnames.get((row.chr1,row.pos1),np.nan),compartmentnames.get((row.chr2,row.pos2),np.nan)
        try:
            mat[int(pos1),int(pos2)] += 1
            mat[int(pos2),int(pos1)] += 1
        except:pass
    if weighted == True :
        complength = comp.counts.value_counts().sort_index().to_numpy()
        mat  = mat/(complength * complength.T)
    return(mat)



def chro_contact_matrix(contacts,chromlist,chrom_len ,resolution = 1000000):
    chro_mat = {}
    contacts = contacts.copy()
    contacts[["pos1","pos2"]] = (contacts[["pos1","pos2"]]/resolution).astype("int")
    chrom_groupby = contacts.groupby(['chr1','chr2'])
    for chrom in chromlist:
        data = chrom_groupby.get_group((chrom,chrom))
        matlen = chrom_len[chrom]
        contact_mat = np.histogram2d(data.pos1,data.pos2,bins=(matlen,matlen))
        chro_mat[chrom] = contact_mat[0]
    return(chro_mat)

