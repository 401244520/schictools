# -*- coding: utf-8 -*-
__author__ = 'zlwang'
import pandas as pd
import gzip,time,os
#from chromatain_demultiplex import chro_contact_matrix
#from matrix_plot import decomposition_mainvector
from .embedding1 import call_compartments,compartment_counts
from .embedding2 import call_compartments2,compartment_concentrate,chro_contact_matrix
## 本模块提供Contact_Matrix这一类，包含交互矩阵，名字，总reads，各染色体的成分等属性
## 并且将三种常见的Hi-C交互数据分别定义为子类进行初始化（格式化）
## Contact_Matrix类目前仅提供了染色体作图和主成分提取的功能，更多拓展功能有待完善


class Contact_Matrix():
    def __init__(self):
        self.contact = None
        self.name = None
        self.amount =None
        self.matrix = None
        self.chrom = {}
        self.resolution = 1000000
        self.compartment = None
        self.compartmentname = None
    
    # 现暂时将降维部分放入main内，本版块仅提供预处理
    def compartment_embedding2(self,resolution = 1000000,refcompartment = "./data/ES-E14_bulk.cis.vecs.tsv_2.csv"):
        print('this singlecell hic file has {} reads'.format(self.amount))
        try :
            compartmentname = pd.read_csv(refcompartment)
        except:
            compartmentname = call_compartments2("/store/zlwang/Workspace/data/4DN/mESC_bulk/mESC_bulk.cis.vecs.tsv")
        time2 = time.time()
        compartments_mat = compartment_concentrate(self.contact,compartmentname,resolution)
        time3 = time.time()
        print('concentrate compa1tments use {} seconds'.format(time3-time2))
        self.compartmentname = compartmentname
        return(compartments_mat)
    # 输出全基因组交互矩阵
    def contact_matrix(self,resolution = 1000000):
        time1 = time.time()
        bins = pd.read_csv("./data/bins_1mb.csv")
        contact_mat = compartment_concentrate(self.contact,bins,resolution)
        print("This singlecell hic file has %r reads \ncalculate the matrix use %.2f seconds (pid %r finished)" % (self.amount,time.time()-time1,os.getpid()))
        self.matrix = contact_mat
        return(contact_mat)
    #
    def plot_matrix(self,chromlist = None,resolution = 1000000):
        chrom = pd.read_csv("./data/chroms.txt",index_col="name")
        chrom = (chrom/resolution).astype(int)
        chrom_len = chrom["length"].to_dict()
        if chromlist == None:
            chromlist = list(chrom_len.keys())
        self.chrom = chro_contact_matrix(self.contact,chromlist = chromlist,chrom_len = chrom_len,resolution = resolution)
        return(self.chrom)
# Type1.eg : chr10	100282258	chr10	132725560
class Contact_type1(Contact_Matrix):
    def __init__(self,path = './data/*type1.txt'):
        Contact_Matrix.__init__(self)
        self.path = path
        self.name = 'type1'
    def preprocess(self):
        try:
            data = pd.read_table(gzip.open(self.path),header=None)
        except:
            data = pd.read_table(self.path,header=None)
        data.columns = ['chr1','pos1','chr2','pos2']
        self.contact = data
        self.amount = len(data)
        print('Type1 data preprocess finished')

# Type2.eg : chr1,chr2,pos1,pos2,strand1,strand2
# 7388231,7388231,3,3,61120032,61120540,False,False
class Contact_type2(Contact_Matrix):
    def __init__(self,path = './data/*type2.csv'):
        Contact_Matrix.__init__(self)
        self.path = path
        self.name = 'type2'
    def preprocess(self):
        try:
            data = pd.read_csv(gzip.open(self.path))
        except:
            data = pd.read_csv(self.path)
        data = data[['chr1','pos1','chr2','pos2']]
        data.chr1 = data.chr1.apply(lambda x : 'chr'+str(x))
        data.chr2 = data.chr2.apply(lambda x : 'chr'+str(x))
        self.contact = data
        self.amount = len(data)
        print('Type2 data preprocess finished')

# Type3.eg : chr1,3007332,1	chr1,3546849,.
class Contact_type3(Contact_Matrix):
    def __init__(self,path = './data/*type3.txt'):
        Contact_Matrix.__init__(self)
        self.path = path
        self.name = 'type3'
    def preprocess(self):
        try:
            data = pd.read_table(gzip.open(self.path),sep ='\t|,',header=None,engine='python')
        except:
            data = pd.read_csv(self.path,sep ='\t|,',header = None,engine='python')
        data.columns = ['chr1','pos1','mark1','chr2','pos2','mark2']
        data = data[['chr1','pos1','chr2','pos2']]
        self.contact = data
        self.amount = len(data)
        print('Type3 data preprocess finished')

