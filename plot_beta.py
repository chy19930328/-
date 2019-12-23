# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 21:25:14 2019

@author: chenhongyu
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.special as ss

def beta_s(x, a, b):
    return x**(a-1)*(1-x)**(b-1)
def beta(x, a, b):
    return beta_s(x, a, b)/ss.beta(a, b)

T_post_alpha=np.array([[45.32,2.02,4.8],
                       [0.39,13.32,3.49],
                       [0.39,3.06,3.2]])
T_post_beta=np.array([[30.21,8.09,19.2],
                      [3.54,5097,18.43],
                      [3.54,517.20,125.48]])
V_post_beta=np.array([[30.21,10.36,14.15],
                       [21.72,99939.9,70.07],
                       [17.18,286.85,64.68]])
V_post_alpha=np.array([[45.32,3.45,3.50],
                       [21.72,4.28,4.61],
                       [11.45,3.2,4.01]])
S_post_alpha=np.array([[71.89,24.01,3.23]])
S_post_beta=np.array([[30.81,158285,26.45]])

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
            
CalculateBetas(Ttrans,Var,Talpha,Tbeta)
CalculateBetas(Vtrans,Var,Valpha,Vbeta)
Calculate_slope_Betas(Slope,Svar,Salpha,Sbeta)

x_01=np.arange(0,1,0.001)
y1=[]
y2=[]
for i in x_01:
    y_1=beta(i,T_post_alpha[2][1],T_post_beta[2][1])
    y1.append(y_1)
    y_2=beta(i,Talpha[2][1],Tbeta[2][1])
    y2.append(y_2)

plt.figure(figsize=(10, 5))
plt.plot(x_01, y1, lw=2)
plt.plot(x_01, y2)

plt.show()



    
