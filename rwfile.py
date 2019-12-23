# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 09:46:49 2019

@author: chenhongyu
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt # plt 用于显示图片
import matplotlib.image as mpimg # mpimg 用于读取图片
import math
 
def loadData(filePath):
    csv = pd.read_csv(filePath, encoding='gbk')
    retData = []
    for line in csv.values:
        retData.append(line[1:].tolist())
        b=np.mat(retData)
    return b

def ThresholdMatrix(vector,f1=0,f2=2):
    row=vector.shape[0]
    column=vector.shape[1]
    for i in range(row):
        for j in range(column):
            if vector[i,j]<f1:
                vector[i,j]=1
            elif vector[i,j]<f2:
                vector[i,j]=2
            else:
                vector[i,j]=3
    return vector

#def MatrixToImage(r,g,b):
    

strImageFileName = "./data/SateliteSmall.jpg"
DataTopology = "./data/RealDataT.csv"    
DataVegetation = "./data/RealDataV.csv"
DataAltitude = "./data/RealDataA.csv"

#显示图片
lena = mpimg.imread(strImageFileName) # 读取和代码处于同一目录下的 lena.png
# 此时 lena 就已经是一个 np.array 了，可以对它进行任意处理
lena.shape #(512, 512, 3)
plt.imshow(lena) # 显示图片
plt.axis('off') # 不显示坐标轴
plt.show()
Tdata=np.flipud(loadData(DataTopology))
Vdata=np.flipud(loadData(DataVegetation))
Adata=np.flipud(loadData(DataAltitude))

mR=ThresholdMatrix(Tdata,50,110)
mG=ThresholdMatrix(Vdata,130,170)
       
class Vertex:
    def __init__(self, key,topology,vegetation,altitude):
        self.id = key
        self.Topology=topology
        self.Vegetation=vegetation
        self.Altitude=altitude
        self.connected_to = {}
        self.probability=0#运行过程中的概率
        self.next_probability=0
        self.total_probability=0#使用转化矩阵生成的固定概率之和
        self.connected_to_p={}#用于计算后验分布使用
        self.local_probability=0
        self.connected_to_next_p={}
    def add_neighbor(self, nbr, weight=0,probability=0):
        self.connected_to[nbr] = weight
        self.connected_to_p[nbr]=probability
        self.connected_to_next_p[nbr]=0
        
    def __str__(self):
        return str(self.id) + ' connected to: ' + str([x.id for x in self.connected_to])
    def get_connections(self):
        return self.connected_to.keys()

    def get_id(self):
        return self.id

    def get_weight(self, nbr):
        return self.connected_to[nbr]

class Graph:
    def __init__(self):
        self.vertex_list = {}
        self.num_vertices = 0

    def add_vertex(self, key,topology,vegetation,altitude):
        self.num_vertices = self.num_vertices + 1
        new_vertex = Vertex(key,topology,vegetation,altitude)
        self.vertex_list[key] = new_vertex
        return new_vertex

    def get_vertex(self, n):
        if n in self.vertex_list:
            return self.vertex_list[n]
        else:
            return None

    def __contains__(self, item):
        return item in self.vertex_list

    def add_edge(self, f, t, cost=0):
        self.vertex_list[f].add_neighbor(t, cost)

    def get_vertices(self):
        return self.vertex_list.keys()

    def __iter__(self):
        return iter(self.vertex_list.values())

def get_elevation(angle):
    if angle>20:
        return 2
    if 0<=angle<20:
        return 1
    if angle<0:
        return 0

def get_mean_T(array):
    A=np.zeros(shape=(34,67))
    array=array[0:170,0:335]
    
    for i in range(0,array.shape[0],5):
        for j in range(0,array.shape[1],5):
            a=np.mean(array[i:i+5,j:j+5])
            x,y=int(i/5),int(j/5)
            A[x,y]=a
    return A
def get_mean_V(array):
    A=np.zeros(shape=(34,67))
    array=array[0:170,0:335]
    for i in range(0,array.shape[0],5):
        for j in range(0,array.shape[1],5):
            section=array[i:i+5,j:j+5].reshape(1,25)
            count_list=np.zeros(3)
            for k in range(25):
                sect=int(section[0,k]-1)
                count_list[sect]+=1
            b=np.argmax(count_list)
            x,y=int(i/5),int(j/5)
            A[x,y]=b+1
    return A

Adata=get_mean_T(Adata)
Vdata=get_mean_V(Vdata)
g = Graph()
row=mR.shape[0]
columns=mR.shape[1]
for i in range(row):
    for j in range(columns):
         g.add_vertex(key=str(i)+'+'+str(j),topology=mR[i,j],vegetation=mG[i,j],altitude=Adata[i,j])

#六边形边长
width=12
#形成包含所有关系的图结构
for i in range(row):
    for j in range(columns):
        if j%2==1:
            if i-1<0:#(i-1,j)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i-1)+'+'+str(j)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i-1)+'+'+str(j),cost=get_elevation(angle*180/math.pi))
            if j-1<0:#(i,j-1)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i)+'+'+str(j-1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i)+'+'+str(j-1),cost=get_elevation(angle*180/math.pi))
            if j+1>=columns:#(i,j+1)有问题
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i)+'+'+str(j+1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i)+'+'+str(j+1),cost=get_elevation(angle*180/math.pi))
            if i+1>=row or j-1<0:#(i+1,j-1)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i+1)+'+'+str(j-1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i+1)+'+'+str(j-1),cost=get_elevation(angle*180/math.pi))
            if i+1>=row:#(i+1,j)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i+1)+'+'+str(j)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i+1)+'+'+str(j),cost=get_elevation(angle*180/math.pi))
            if i+1>=row or j+1>=columns:#(i+1,j+1)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i+1)+'+'+str(j+1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i+1)+'+'+str(j+1),cost=get_elevation(angle*180/math.pi)) 
        if j%2==0:
            if i-1<0:#(i-1,j)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i-1)+'+'+str(j)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i-1)+'+'+str(j),cost=get_elevation(angle*180/math.pi))
            if i-1<0 or j+1>=columns:#(i-1,j+1)
                 pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i-1)+'+'+str(j+1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i-1)+'+'+str(j+1),cost=get_elevation(angle*180/math.pi))
            if  j+1>=columns:#(i,j+1)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i)+'+'+str(j+1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i)+'+'+str(j+1),cost=get_elevation(angle*180/math.pi))
            if  i+1>=row:#(i+1,j)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i+1)+'+'+str(j)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i+1)+'+'+str(j),cost=get_elevation(angle*180/math.pi))
            if  j-1<0:#(i,j-1)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i)+'+'+str(j-1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i)+'+'+str(j-1),cost=get_elevation(angle*180/math.pi))
            if  i-1<0 or j-1<0:#(i-1,j-1)
                pass
            else:
                angle=math.atan((g.vertex_list[str(i)+'+'+str(j)].Altitude-g.vertex_list[str(i-1)+'+'+str(j-1)].Altitude)/(width*math.sqrt(3)))
                g.add_edge(str(i)+'+'+str(j),str(i-1)+'+'+str(j-1),cost=get_elevation(angle*180/math.pi))

np.save('Vdata.npy',Vdata)
print(g.vertex_list['19+16'].connected_to)            
            
          
       
  



            