# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 09:40:01 2019

@author: chenhongyu
"""

import numpy as np
import rwfile
from scipy.stats import beta
import re
import math

Ttrans=np.array([[0.6,0.2,0.2],
                 [0.1,0.7,0.2],
                 [0.1,0.5,0.4]])

Vtrans=np.array([[0.6,0.25,0.15],
                [0.5,0.3,0.2],
                [0.4,0.4,0.2]])  
         
Slope=np.array([0.7,0.2,0.1])

Var=np.array([[0.14*0.14,0.15*0.15,0.1*0.1],
             [0.15*0.15,0.15*0.15,0.15*0.15],
             [0.15*0.15,0.15*0.15,0.15*0.15]])
Svar=np.full((1,3),0.15*0.15)

Tbeta,Talpha,Vbeta,Valpha=np.zeros(shape=(3,3)),np.zeros(shape=(3,3)),np.zeros(shape=(3,3)),np.zeros(shape=(3,3))
Sbeta,Salpha=np.zeros(shape=(1,3)),np.zeros(shape=(1,3))

def CalculateBetas(means,varition,alphas,betas):
    row=means.shape[0]
    column=means.shape[1]
    for i in range(row):
        for j in range(column):
            mean = means[i, j]
            var = varition[i,j]
            temp = (1 - mean) / mean
            alpha = (1 - var * temp) / (var * temp * (1 + temp))
            beta = alpha * temp
            alphas[i, j] = alpha;
            betas[i, j] = beta;
            
def Calculate_slope_Betas(means,varition,alphas,betas):
    row=means.shape[0]
    for i in range(row):
            mean = means[i]
            var = varition[0,i]
            temp = (1 - mean) / mean
            alpha = (1 - var * temp) / (var * temp * (1 + temp))
            beta = alpha * temp
            alphas[0,i] = alpha;
            betas[0,i] = beta;


def sample_prior(alpha,beta):
    A=np.zeros(shape=(3,3))
    for i in range(3):
        for j in range(3):
            A[i,j]=np.random.beta(alpha[i,j],beta[i,j])
    for i in range(3):
        for j in range(3):
            A[i,j]=A[i,j]/sum(A[i])
    return A

def sample_prior_slope(alpha,beta):
    A=np.zeros(shape=(1,3))
    for i in range(3):
        A[0,i]=np.random.beta(alpha[0,i],beta[0,i])
    for i in range(3):
        A[0,i]=A[0,i]/sum(A[0])        
    return A

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
    states =[cur]
    for i in range(1000000):
        nex, u = np.random.rand(),np.random.rand()
        if u < np.min((normal_s(nex, a, b)/normal_s(cur, a, b), 1)):
            states.append(nex)
            if len(states)==10000:
                break
                cur = nex    
    return states
def Calculate_post_Betas(array):
    mean = np.mean(array)
    var = np.var(array)
    temp = (1 - mean) / mean
    alpha = (1 - var * temp) / (var * temp * (1 + temp))
    beta = alpha * temp
    return alpha,beta
              
#生成三个转换矩阵的alpha、beta分布
CalculateBetas(Ttrans,Var,Talpha,Tbeta)
CalculateBetas(Vtrans,Var,Valpha,Vbeta)
Calculate_slope_Betas(Slope,Svar,Salpha,Sbeta)
g=rwfile.g
data_raw=loadData('data.old.txt',0.8)

for k in range(2):
    for  i in data_raw:
        Ttrans_sample=sample_prior(Talpha,Tbeta)
        Vtrans_sample=sample_prior(Valpha,Vbeta)
        Slope_sample=sample_prior_slope(Salpha,Sbeta)
        location_in=i[0]+'+'+i[1]
        location_out=i[2]+'+'+i[3]
        Trow_update=int(g.vertex_list[location_in].Topology)-1
        Tcolumn_update=int(g.vertex_list[location_out].Topology)-1
        Vrow_update=int(g.vertex_list[location_in].Vegetation)-1
        Vcolumn_update=int(g.vertex_list[location_out].Vegetation)-1
        slope_update=int(g.vertex_list[location_in].connected_to[location_out])
        #得到三个参数的10000个值
        T_sample=mcmc_normal(Ttrans[Trow_update,Tcolumn_update],Var[Trow_update,Tcolumn_update])
        V_sample=mcmc_normal(Vtrans[Vrow_update,Vcolumn_update],Var[Vrow_update,Vcolumn_update])
        S_sample=mcmc_normal(Slope[slope_update],Svar[0,slope_update])

        #得到一个点的各个方向的概率
        for v in g.vertex_list[location_in].connected_to.items():
            Trow=int(g.vertex_list[location_in].Topology)-1
            Tcolumn=int(g.vertex_list[v[0]].Topology)-1
            Vrow=int(g.vertex_list[location_in].Vegetation)-1
            Vcolumn=int(g.vertex_list[v[0]].Vegetation)-1
            slope=int(g.vertex_list[location_in].connected_to[v[0]])
            g.vertex_list[location_in].connected_to_p[v[0]]=Ttrans_sample[Trow,Tcolumn]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]
            g.vertex_list[location_in].local_probability=Ttrans_sample[Trow,Trow]*Ttrans_sample[Vrow,Vrow]*Slope_sample[0,1]
            g.vertex_list[location_in].total_probability+=g.vertex_list[location_in].connected_to_p[v[0]]
        g.vertex_list[location_in].total_probability+=g.vertex_list[location_in].local_probability
        reduce_value=g.vertex_list[location_in].connected_to_p[location_out]
        T_likelihood_list=[]
        V_likelihood_list=[]
        S_likelihood_list=[]
        #计算一万个取值的总概率变化列表
        #更新T
        for i in range(10000):
            Ttrans_sample[Trow_update][Tcolumn_update]=T_sample[i]
            Ttrans_sample[Trow_update]=Ttrans_sample[Trow_update]/sum(Ttrans_sample[Trow_update])
            totalprobability=g.vertex_list[location_in].total_probability-reduce_value+T_sample[i]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]
            g.vertex_list[location_in].connected_to_p[location_out]=T_sample[i]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]
            local_probability=g.vertex_list[location_in].local_probability/totalprobability
            T_likelihood_list.append(g.vertex_list[location_in].connected_to_p[location_out])
        Talpha[Trow_update][Tcolumn_update],Tbeta[Trow_update][Tcolumn_update]=Calculate_post_Betas(T_likelihood_list)    
        print('-------------')
        ##更新V
        for i in range(10000):         
            totalprobability=g.vertex_list[location_in].total_probability-reduce_value+Ttrans_sample[Vrow,Vcolumn]*V_sample[i]*Slope_sample[0,slope]
            g.vertex_list[location_in].connected_to_p[location_out]=Ttrans_sample[Vrow,Vcolumn]*V_sample[i]*Slope_sample[0,slope]
            local_probability=g.vertex_list[location_in].local_probability/totalprobability
            V_likelihood_list.append(g.vertex_list[location_in].connected_to_p[location_out])
        Valpha[Trow_update][Tcolumn_update],Vbeta[Trow_update][Tcolumn_update]=Calculate_post_Betas(V_likelihood_list)
        ##更新s
        for i in range(10000):
            totalprobability=g.vertex_list[location_in].total_probability-reduce_value+Ttrans_sample[Vrow,Vcolumn]*Vtrans_sample[Vrow,Vcolumn]*S_sample[i]
            g.vertex_list[location_in].connected_to_p[location_out]=Ttrans_sample[Vrow,Vcolumn]*Vtrans_sample[Vrow,Vcolumn]*S_sample[i]
            local_probability=g.vertex_list[location_in].local_probability/totalprobability
            S_likelihood_list.append(g.vertex_list[location_in].connected_to_p[location_out])
        Salpha[0][Tcolumn_update],Sbeta[0][Tcolumn_update]=Calculate_post_Betas(S_likelihood_list)
    









