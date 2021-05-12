# -*- coding: utf-8 -*-
# 本模块提供染色体的交互矩阵的嵌入表示
# 目前只针对一个染色体进行计算,计算染色体上的覆盖度

import pandas as pd
import numpy as np
import random

def call_compartments(compartmentpath,resolution = 1000000):
    # TODO：use cooltools api rather than cli to calculate comparment engivectors
    data = pd.read_table(compartmentpath)
    c = ['chrom',0,0,1,0]
    new = []
    for i in data.index:
        d = data.loc[i].copy()
        if d[4] == np.nan:  #使用E2作为划分标准
            c[1] = d[1]
        elif c[4]*d[4]>= 0 :
            c[2] = d[2]
            c[4] += d[4]
        elif d[4]*c[4]<0 :
            new.append(list(c))
            c = d       
        else :
            pass
    compartments = pd.DataFrame(new)
    compartments = compartments.drop(index = 0) #是否去除染色体末端空值
    compartments.columns = ['chrom','start','end','vector','length','counts']
    compartments[['start','end']] = (compartments[['start','end']]/resolution).astype("int")
    compartments['length'] = compartments.end - compartments.start
    compartments["counts"] = 0
    compartments.to_csv("{}.csv".format(compartmentpath),index=0)
    return(compartments)

def compartment_counts(contacts,compartments,contact_store = 'all',resolution = 1000000):
    print(contact_store)
    if contact_store == "all":
        pass
    elif contact_store == "intra":
        contacts = contacts[contacts.chr1 == contacts.chr2]
    elif contact_store == "inter":
        contacts = contacts[contacts.chr1 != contacts.chr2]
    contact1 =contacts[["chr1","pos1"]]
    contact1.columns = ["chr","pos"]
    contact2 =contacts[["chr2","pos2"]]
    contact2.columns = ["chr","pos"]
    newcontacts = pd.concat([contact1,contact2],ignore_index=True)
    newcontacts["pos"] = (newcontacts["pos"]/resolution).astype("int")
    group = newcontacts.groupby("chr")
    groupcounts = group["pos"].value_counts()
    embedding = []
    for ind,data in compartments.iterrows():
        counts = 0
        for i in range(data.start,data.end+1): #双边界
            try :
                counts += groupcounts[data.chrom,i]
            except:
                pass
        embedding.append(counts)
    return(embedding)


