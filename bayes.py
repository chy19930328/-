# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:44:51 2019

@author: chenhongyu
"""
import numpy as np
import rwfile
from scipy.stats import beta
from matplotlib import pyplot as plot #用来绘制图形
from mpl_toolkits.mplot3d import Axes3D  #用来给出三维坐标系。
#根据专家意见得到的先验
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
               
#生成三个转换矩阵的alpha、beta分布
CalculateBetas(Ttrans,Var,Talpha,Tbeta)
CalculateBetas(Vtrans,Var,Valpha,Vbeta)
Calculate_slope_Betas(Slope,Svar,Salpha,Sbeta)
#class ExpertOpinions:
#    def __init__(self,beta_topology,alpha_vegetable,beta_vegetation,beta_vegetation,slope):
#        self.alphasT=np.zeros(shape=(3,3))
#        self.betasT=np.zeros(shape=(3,3))
#        self.alphasV=np.zeros(shape=(3,3))
#        self.betasV=np.zeros(shape=(3,3))
#        self.alpasS=np.zeros(shape=(1,3))
#        self.betasS=np.zeros(shape=(1,3))
#        for i in range(topology.shape[0]):
#            for j in range(topology.shape[1]):
#                pass

g=rwfile.g
#计算每个点的走出概率总和，用于概率归一化

def vertex_compute_probability(Ttrans_sample,Vtrans_sample,Slope_sample):
    for k in g.vertex_list.items():
        for v in g.vertex_list[k[0]].connected_to.items():
            Trow=int(g.vertex_list[k[0]].Topology)-1
            Tcolumn=int(g.vertex_list[v[0]].Topology)-1
            Vrow=int(g.vertex_list[k[0]].Vegetation)-1
            Vcolumn=int(g.vertex_list[v[0]].Vegetation)-1
            slope=int(g.vertex_list[k[0]].connected_to[v[0]])
            P_choose=Ttrans_sample[Trow,Tcolumn]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]
            g.vertex_list[k[0]].total_probability+=P_choose
        P_local_choose=Ttrans_sample[Trow,Trow]*Ttrans_sample[Vrow,Vrow]*Slope_sample[0,1]
        g.vertex_list[k[0]].total_probability+=P_local_choose
    for k in g.vertex_list.items():
        for v in g.vertex_list[k[0]].connected_to.items():
            Trow=int(g.vertex_list[k[0]].Topology)-1
            Tcolumn=int(g.vertex_list[v[0]].Topology)-1
            Vrow=int(g.vertex_list[k[0]].Vegetation)-1
            Vcolumn=int(g.vertex_list[v[0]].Vegetation)-1
            slope=int(g.vertex_list[k[0]].connected_to[v[0]])
            g.vertex_list[k[0]].connected_to_p[v[0]]=Ttrans_sample[Trow,Tcolumn]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]/g.vertex_list[k[0]].total_probability
        g.vertex_list[k[0]].local_probability=Ttrans_sample[Trow,Trow]*Ttrans_sample[Vrow,Vrow]*Slope_sample[0,1]/g.vertex_list[k[0]].total_probability
        

location=input('请输入最后发现失踪人员位置信息：')
Time=int(input('请输入已经失踪的时长：'))
#def bulid_path(graph,location,time):
#    pass
t=0
Ploc=1#最初位置的概率为1
Path=[]
Path.append(location)
while t<=Time:
    Ttrans_sample=sample_prior(Talpha,Tbeta)
    Vtrans_sample=sample_prior(Valpha,Vbeta)
    Slope_sample=sample_prior_slope(Salpha,Sbeta)
    path_probability=[]
    Pmax=0
    total_p=0
        #从beta分布中更新采样
    Ttrans_sample=sample_prior(Talpha,Tbeta)
    Vtrans_sample=sample_prior(Valpha,Vbeta)
    Slope_sample=sample_prior_slope(Salpha,Sbeta)
    for kv in g.vertex_list[location].connected_to.items():
        Trow=int(g.vertex_list[location].Topology)-1
        Tcolumn=int(g.vertex_list[kv[0]].Topology)-1
        Vrow=int(g.vertex_list[location].Vegetation)-1
        Vcolumn=int(g.vertex_list[kv[0]].Vegetation)-1
        slope=int(g.vertex_list[location].connected_to[kv[0]])
        P_choose=Ploc*Ttrans_sample[Trow,Tcolumn]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]
        P_local_choose=Ploc*Ttrans_sample[Trow,Trow]*Ttrans_sample[Vrow,Vrow]*Slope_sample[0,1]
        
        #计算出那个临近点为概率最大
        if Pmax<P_choose:
            Pmax=P_choose
            next_part=kv[0]
        path_probability.append(P_choose)
    if Pmax<P_local_choose:
        next_part=location
    path_probability.append(P_local_choose)
    path_probability=[i/sum(path_probability) for i in path_probability]
    Path.append(location)
    Path.append(next_part)
    location=next_part
    t+=1
    
#print(Path)
def build_probabilitymatrixmap(location,Time):
    Ttrans_sample=sample_prior(Talpha,Tbeta)
    Vtrans_sample=sample_prior(Valpha,Vbeta)
    Slope_sample=sample_prior_slope(Salpha,Sbeta)
    vertex_compute_probability(Ttrans_sample,Vtrans_sample,Slope_sample)  
    t=0
    g.vertex_list[location].probability=1
    while t<=Time:
        for k in g.vertex_list.keys():#遍历图中所有点
            g.vertex_list[k].next_probability=0
            #遍历某一点的所有邻接点，计算邻接点到这一点的概率和，就是下一秒的到达这一点概率
            for kv in g.vertex_list[k].connected_to.items():
                Trow=int(g.vertex_list[kv[0]].Topology)-1
                Tcolumn=int(g.vertex_list[k].Topology)-1
                Vrow=int(g.vertex_list[kv[0]].Vegetation)-1
                Vcolumn=int(g.vertex_list[k].Vegetation)-1
                slope=int(g.vertex_list[kv[0]].connected_to[k])
                if g.vertex_list[kv[0]].probability>0:
                #计算邻接点到该点的概率
                    g.vertex_list[k].next_probability+=g.vertex_list[kv[0]].probability*Ttrans_sample[Trow,Tcolumn]*Vtrans_sample[Vrow,Vcolumn]*Slope_sample[0,slope]/g.vertex_list[kv[0]].total_probability
            #计算留在该点的概率
            P_local_choose=g.vertex_list[k].probability*Ttrans_sample[Tcolumn,Tcolumn]*Ttrans_sample[Vcolumn,Vcolumn]*Slope_sample[0,1]/g.vertex_list[k].total_probability
            g.vertex_list[k].next_probability+=P_local_choose
                #最终得到这一轮生成的概率
                #更新概率，为下一轮迭代做准备        
        for k in g.vertex_list.keys():
            g.vertex_list[k].probability=g.vertex_list[k].next_probability            
        t+=1

build_probabilitymatrixmap(location,Time)
figure = plot.figure() 
#画出三维坐标系：
axes = Axes3D(figure)
X = np.arange(0, 34, 1)
Y = np.arange(0, 67, 1)
#限定图形的样式是网格线的样式：
Z=np.zeros(shape=(67,34))
for i in range(67):
    for j in range(33):
        Z[i,j]=g.vertex_list[str(j)+'+'+str(i)].probability
X, Y = np.meshgrid(X, Y)#绘制曲面，采用彩虹色着色：
axes.plot_surface(X, Y, Z,cmap='rainbow')
#图形可视化：
plot.show()