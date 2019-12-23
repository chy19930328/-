# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 12:22:49 2019

@author: chenhongyu
"""
import numpy as np
import rwfile
import bayes
from matplotlib import pyplot as plot #用来绘制图形
from mpl_toolkits.mplot3d import Axes3D  #用来给出三维坐标系。

g=rwfile.g
Ttrans_sample=bayes.Ttrans_sample
Vtrans_sample=bayes.Vtrans_sample
Slope_sample=bayes.Slope_sample

location=bayes.location
Time=bayes.Time

#计算生成概率地图矩阵
def build_probabilitymatrixmap(location,Time):
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