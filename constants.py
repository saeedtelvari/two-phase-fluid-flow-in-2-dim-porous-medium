#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


properties = pd.read_excel('Property.xlsx')


# In[3]:


properties


# In[12]:


dx = np.array(pd.read_excel('dx.xlsx', header = None))
dy = np.array(pd.read_excel('dy.xlsx', header = None))
dz = np.array(pd.read_excel('dz.xlsx', header = None))


# In[5]:


BCW = np.array(pd.read_excel('BCW.xlsx', header = None))
BCN = np.array(pd.read_excel('BCN.xlsx', header = None))
BCS = np.array(pd.read_excel('BCS.xlsx', header = None))
BCE = np.array(pd.read_excel('BCE.xlsx', header = None))
IB = properties['IBC(1=Constant Flow Rate 2=Constant BHP)'][0]


# In[3]:


ndt = 10000
# Wells
BHPo = np.zeros([ndt+1, ndx, ndy])
BHPw = np.zeros([ndt+1, ndx, ndy])
qo = np.zeros([ndt+1, ndx, ndy])
qw = np.zeros([ndt+1, ndx, ndy])
rw_o = np.zeros([ndt+1, ndx, ndy])
rw_w = np.zeros([ndt+1, ndx, ndy])
So = np.zeros([ndt+1, ndx, ndy])
Sw = np.zeros([ndt+1, ndx, ndy])
BHPo[0] = np.array(pd.read_excel('BHPo.xlsx', header = None))
BHPw[0] = np.array(pd.read_excel('BHPw.xlsx', header = None))
qo[0] = np.array(pd.read_excel('qo.xlsx', header = None))
qw[0] = np.array(pd.read_excel('qw.xlsx', header = None))
rw_o[0] = np.array(pd.read_excel('rwo.xlsx', header = None))
rw_w[0] = np.array(pd.read_excel('rww.xlsx', header = None))
So[0] = np.array(pd.read_excel('SKw.xlsx', header = None))
Sw[0] = np.array(pd.read_excel('SKo.xlsx', header = None))
Bo_0 = properties['Bo0'][0]
Bw_0 = properties['Bw0'][0]
Ro_o = properties['ρo'][0]
Ro_w = properties['ρw'][0]
#
dt = properties['dt'][0]
ndt = properties['ndt'][0]
Poi = properties['Poi'][0]
Co = properties['co'][0]
Cw = properties['cw'][0]
Cf = properties['cf'][0]
P_0 = properties['p0'][0]
IBo = properties.iloc[0][12]
IBw = properties.iloc[0][13]

Kx = np.array(pd.read_excel('Kx.xlsx', header = None))
Ky = np.array(pd.read_excel('Ky.xlsx', header = None))
porosity_0 = pd.read_excel('phi0.xlsx', header = None)
rel_ro_o = 0.00021584 * Ro_o * 32.17
rel_ro_w = 0.00021584 * Ro_w * 32.17
ndx, ndy = dx.shape

r_eq = 0.14 * np.power((np.power(dx, 2) + np.power(dy, 2)), 0.5) # Equivalent drainage radius for each grid block 


# In[10]:


# Calclates Area and then the Transmissbility of gridblocks
def H(BCW, BCN, BCE, BCS, dx, dy, dz, Kx, Ky, ndx, ndy):
    Ax = np.zeros([ndx, ndy])
    Ay = np.zeros([ndx, ndy])
    h = np.zeros([4, ndx, ndy])
    for j in range(ndy):
        for i in range(ndy):
            Ax[i][j] = dx[i][j] * dz[i][j]
            Ay[i][j] = dy[i][j] * dz[i][j]
    for j in range(ndy):
        for i in range(ndx):
            # North
            if BCN[i][j] == 100:
                h[0][i][j] = 2 * Ax[i][j] * Ax[i-1][j] * Kx[i][j] * Kx[i-1][j] /                 (Ax[i][j] * Kx[i][j] * dx[i-1][j] + Ax[i-1][j] * Kx[i-1][j] * dx[i][j]) 
            else:
                h[0][i][j] = Ax[i][j] * Kx[i][j] / dx[i][j]
            # West
            if BCW[i][j] == 100:
                h[1][i][j] = 2 * Ay[i][j] * Ay[i][j-1] * Ky[i][j] * Ky[i][j-1] /                 (Ay[i][j] * Ky[i][j] * dy[i][j-1] + Ay[i][j-1] * Ky[i][j-1] * dy[i][j]) 
            else:
                h[1][i][j] = Ay[i][j] * Ky[i][j] / dy[i][j]
            # East
            if BCE[i][j] == 100:
                h[2][i][j] = 2 * Ay[i][j] * Ay[i][j+1] * Ky[i][j] * Ky[i][j+1] /                 (Ay[i][j] * Ky[i][j] * dy[i][j+1] + Ay[i][j+1] * Ky[i][j+1] * dy[i][j]) 
            else:
                h[2][i][j] = Ay[i][j] * Ky[i][j] / dy[i][j]
            # South
            if BCS[i][j] == 100:
                h[3][i][j] = 2 * Ax[i][j] * Ax[i+1][j] * Kx[i][j] * Kx[i+1][j] /                 (Ax[i][j] * Kx[i][j] * dx[i+1][j] + Ax[i+1][j] * Kx[i+1][j] * dx[i][j]) 
            else:
                h[3][i][j] = Ax[i][j] * Kx[i][j] / dx[i][j]
    return h[0], h[1], h[2], h[3]


# In[11]:


HN, HW, HE, HS= H(BCW, BCN, BCE, BCS, dx, dy, dz, Kx, Ky, ndx, ndy)

