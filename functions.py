#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


# Functions :
# relative_permeability
# mobility
# viscosity
# FVF
# capillary_pressure


# In[7]:


def relative_permeability(S_wr, S_or, S_w, ndx, ndy):
    S_we = np.zeros([ndx, ndy])
    Kr_o = np.zeros([ndx, ndy])
    Kr_w = np.zeros([ndx, ndy])
    
#     Kr_o = ((-S_w + 1 - S_or)/(1 - S_or - S_wr)) ** 6
#     Kr_w = ((S_w - S_wr)/(1 - S_or - S_wr)) ** 6
    S_we = (np.subtract(S_w, S_wr)) / (1 - S_or - S_wr)
    Kr_o = 0.9 * np.power(1-S_we, 2)
    Kr_w = 0.6 * np.power(S_we, 2)   
    for i in range(ndx):
        for j in range(ndy):
            if Kr_o[i][j] > 1:
                Kr_o[i][j] = 1
            if Kr_w[i][j] > 1:
                Kr_w[i][j] = 1

    return Kr_o, Kr_w


# In[3]:


def mobility(BCW, BCE, BCN, BCS, Kr_o, Kr_w, Vis_o, Vis_w, Bo, Bw, Po, Pw, ndx, ndy):
    lambda_o = np.zeros([5, ndx, ndy])
    lambda_w = np.zeros([5, ndx, ndy])
    for i in range(ndx):
        for j in range(ndy):
            # West
            if BCW[i][j] == 100:
                if Po[i][j] > Po[i][j-1]:
                    lambda_o[1][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                else:
                    lambda_o[1][i][j] = Kr_o[i][j-1] / (Vis_o[i][j-1] * Bo[i][j-1])
                if Pw[i][j] > Pw[i][j-1]:
                    lambda_w[1][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
                else:
                    lambda_w[1][i][j] = Kr_w[i][j-1] / (Vis_w[i][j-1] * Bw[i][j-1])
            else:
                lambda_o[1][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                lambda_w[1][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
            # East
            if BCE[i][j] == 100:
                if Po[i][j] > Po[i][j+1]:
                    lambda_o[2][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                else:
                    lambda_o[2][i][j] = Kr_o[i][j+1] / (Vis_o[i][j+1] * Bo[i][j+1])
                if Pw[i][j] > Pw[i][j+1]:
                    lambda_w[2][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
                else:
                    lambda_w[2][i][j] = Kr_w[i][j+1] / (Vis_w[i][j+1] * Bw[i][j+1])
            else:
                lambda_o[2][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                lambda_w[2][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
            # North
            if BCN[i][j] == 100:
                if Po[i][j] > Po[i-1][j]:
                    lambda_o[3][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                else:
                    lambda_o[3][i][j] = Kr_o[i-1][j] / (Vis_o[i-1][j] * Bo[i-1][j])
                if Pw[i][j] > Pw[i-1][j]:
                    lambda_w[3][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
                else:
                    lambda_w[3][i][j] = Kr_w[i-1][j] / (Vis_w[i-1][j] * Bw[i-1][j])
            else:
                lambda_o[3][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                lambda_w[3][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
            # South
            if BCS[i][j] == 100:
                if Po[i][j] > Po[i+1][j]:
                    lambda_o[4][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                else:
                    lambda_o[4][i][j] = Kr_o[i+1][j] / (Vis_o[i+1][j] * Bo[i+1][j])
                if Pw[i][j] > Pw[i+1][j]:
                    lambda_w[4][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
                else:
                    lambda_w[4][i][j] = Kr_w[i+1][j] / (Vis_w[i+1][j] * Bw[i+1][j])
            else:
                lambda_o[4][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j])
                lambda_w[4][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j])
            
            lambda_o[0][i][j] = Kr_o[i][j] / (Vis_o[i][j] * Bo[i][j]) # for first initialization of lambda
            lambda_w[0][i][j] = Kr_w[i][j] / (Vis_w[i][j] * Bw[i][j]) # for first initialization of lambda
    return lambda_o, lambda_w


# In[2]:


def viscosity(Po, Pw, ndx, ndy):
    # We could use pressure dependent equations but I just decided to use constant variable for each of the viscosities
    # I will add the equations later
    Vis_o = np.zeros([ndx, ndy])
    Vis_w = np.zeros([ndx, ndy])
    Vis_o += 2
    Vis_w += 2
    return Vis_o, Vis_w


# In[7]:


def FVF(Po, Pw, Bo_0, Bw_0, Co, Cw, P_0, ndx, ndy, t):
    Bo = np.zeros([ndx, ndy])
    Bw = np.zeros([ndx, ndy])
    dBo = np.zeros([ndx, ndy])
    dBw = np.zeros([ndx, ndy])
    Bo = Bo_0 / (1 + Co * (Po - P_0))
    Bw = Bw_0 / (1 + Cw * (Pw - P_0))
    for i in range(ndx):
        for j in range(ndy):
            dBo[i][j] = Co / Bo_0
            dBw[i][j] = Cw / Bw_0
    return Bo, Bw, dBo, dBw


# In[5]:


def capillary_pressure(S_wr, S_or, S_w, ndx, ndy):
    S_we = np.zeros([ndx, ndy])
    Pc = np.zeros([ndx, ndy]) 
    dPc = np.zeros([ndx, ndy]) 
    
    
    for i in range(ndx):
        for j in range(ndy):
            S_we[i][j] = (S_w[i][j] - S_wr) / (1 - S_wr - S_or)
            if S_we[i][j] < 0.01:
                Pc[i][j] = 150
                dPc[i][j] = -15/3 * (1/(1 - S_wr - S_or)) / (np.power(0.001, 4/3))
            else:
                Pc[i][j] = 15 / np.power(S_we[i][j], 1/3)
                dPc[i][j] = -15/3 * (1/(1 - S_wr - S_or)) / (np.power(S_we[i][j], 4/3))
    return Pc, dPc


# In[8]:


def porosity(P, Cphi, phi_0, P_0, ndx, ndy):
    phi = np.zeros([ndx, ndy])
    for i in range(ndx):
        for j in range(ndy):
            phi[i][j] = phi_0[i][j] * (1 + Cphi * (P[i][j] - P_0))
    return phi

