# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 20:19:08 2019

@author: chenhongyu
"""
import numpy as np
import scipy.special as ss
import bayes
import math
import rwfile

Talpha,Tbeta,Valpha,Vbeta,Salpha,Sbeta=bayes.Talpha,bayes.Tbeta,bayes.Valpha,bayes.Vbeta,bayes.Salpha,bayes.Sbeta
g=rwfile.g   
#使用自己的方法更新beta分布概率
def update_Ttrans(means,alphas,betas,loaction,nextpart):
    row=int(g.vertex_list[location].Topology)-1
    column=int(g.vertex_list[nextpart].Topology)-1
    alphas[row,column]+=1
    for i in range(betas.shape[0]):
        for j in range(alphas.shape[1]):
            betas[i,j]+=1
            #Mean(‘m’), variance(‘v’), skew(‘s’), and/or kurtosis(‘k’).
            mean, var, skew, kurt = beta.stats(alphas[i,j], betas[i,j], moments='mvsk')
            means[i,j]=mean
def update_Vtrans(means,alphas,betas,loaction,nextpart):
    row=int(g.vertex_list[location].Vegetation)-1
    column=int(g.vertex_list[nextpart].Vegetation)-1
    alphas[row,column]+=1
    for i in range(betas.shape[0]):
        for j in range(alphas.shape[1]):
            betas[i,j]+=1
            #Mean(‘m’), variance(‘v’), skew(‘s’), and/or kurtosis(‘k’).
            mean, var, skew, kurt = beta.stats(alphas[i,j], betas[i,j], moments='mvsk')
            means[i,j]=mean

def loadData(filename,split,trainingSet=[],testSet=[]):
    with open(filename,'r') as file:
        lines = file.readlines()
        dataset = [[] for i in range(len(lines)-1)]
        for i in range(len(dataset)):
            dataset[i][:] = (item for item in lines[i].strip().split(','))   # 逐行读取数据        print("dateset:",dataset,len(dataset))
    return dataset
    
def beta_s(x, a, b):
    return x**(a-1)*(1-x)**(b-1)
def beta(x, a, b):
    return beta_s(x, a, b)/ss.beta(a, b)
def normal_s(x,u,s):
    sig=math.sqrt(s)
    return np.exp(-(x - u) ** 2 /(2* sig **2))/(math.sqrt(2*math.pi)*sig)

def mcmc_beta(a, b):
    cur = np.random.rand()
    states = [cur]
    for i in range(100000):
        nex, u = np.random.rand(),np.random.rand()
        if u < np.min((beta_s(nex, a, b)/beta_s(cur, a, b), 1)):
            states.append(nex)
            if len(states)>=10000:
                break
            cur = nex
    mean=np.mean(states)
    var=np.var(states)
    return mean,var

def mcmc_normal(a, b):
    cur = np.random.rand()
    states = [cur]
    for i in range(100000):
        nex, u = np.random.rand(),np.random.rand()
        if u < np.min((normal_s(nex, a, b)/normal_s(cur, a, b), 1)):
            states.append(nex)
            if len(states)>=10000:
                break
            cur = nex
    mean=np.mean(states)
    var=np.var(states)
    
    return states

data_raw=loadData('data.old.txt',0.8)
post_data_location=[[i[0]+'+'+i[1],i[2]+'+'+i[3]] for i in data_raw]#点对表示从这里到那里
T_post_location=np.zeros(shape=(134,2))
k=0
for i in post_data_location:
    T_post_location[k,0]=g.vertex_list[i[0]].Topology
    T_post_location[k,1]=g.vertex_list[i[1]].Topology
    k+=1
T_sample=[]
V_sample=[]
for _ in range(1):
    alpha=Talpha[2,2]
    beta=Tbeta[2,2]
    Vrow=Valpha[2,2]
    Vcolumn=Vbeta[2,2]
    T_sample.append(mcmc_normal(alpha,beta))
    V_sample.append(mcmc_normal(Vrow,Vcolumn))
    
#def calcalate_post(alphas,betas):
#    means=np.zeros(shape=(3,3))
#    var=means=np.zeros(shape=(3,3))
#    for i in range(3):
#        for j in range(3):
#            means[i,j],var[i,j]=mcmc_beta(alphas[i,j],betas[i,j])
#            for _ in range(1):
#                means[i,j],var[i,j]=mcmc_normal(means[i,j], var[i,j])
#    return means,var
#
#def calcalate_post_s(alphas,betas):
#    means=np.zeros(shape=(1,3))
#    var=means=np.zeros(shape=(1,3))
#    for j in range(3):
#        means[0,j],var[0,j]=mcmc_beta(alphas[0,j],betas[0,j])
#        for _ in range(1):
#            means[0,j],var[0,j]=mcmc_normal(means[0,j], var[0,j])
#    return means,var  
def sample_prior(alpha,beta):
    A=np.zeros(shape=(3,3,2))
    for i in range(3):
        for j in range(3):
            A[i,j]=mcmc_beta(alpha[i,j],beta[i,j])
    return A

def sample_prior_slope(alpha,beta):
    A=np.zeros(shape=(1,3,2))
    for i in range(3):
        A[0,i]=mcmc_beta(alpha[0,i],beta[0,i])
    return A
#使用吉布斯采样得到beta分布的先验
Ttrans_sample=sample_prior(Talpha,Tbeta)
Vtrans_sample=sample_prior(Valpha,Vbeta)
Slope_sample=sample_prior_slope(Salpha,Sbeta)                    
location_now,location_future='22+9','21+9'
path_probability=[]
Ploc=1
Path=[]
for kv in g.vertex_list[location_now].connected_to.items():
    Trow=int(g.vertex_list[location_now].Topology)-1
    Tcolumn=int(g.vertex_list[kv[0]].Topology)-1
    Vrow=int(g.vertex_list[location_now].Vegetation)-1
    Vcolumn=int(g.vertex_list[kv[0]].Vegetation)-1
    slope=int(g.vertex_list[location_now].connected_to[kv[0]])
    P_choose=Ploc*Ttrans_sample[Trow,Tcolumn][0]*Vtrans_sample[Vrow,Vcolumn][0]*Slope_sample[0,slope][0]
    P_local_choose=Ploc*Ttrans_sample[Trow,Trow][0]*Ttrans_sample[Vrow,Vrow][0]*Slope_sample[0,1][0]
    path_probability.append(P_choose)
    Path.append(kv[0])
        #计算出那个临近点为概率最大              
path_probability.append(P_local_choose) 
Path.append(location_now)
#归一化后的7个概率
path_probability_normal=[i/np.sum(path_probability)for i in path_probability]    
path_dict=dict(zip(Path,path_probability_normal))













    
            
    
        
        
