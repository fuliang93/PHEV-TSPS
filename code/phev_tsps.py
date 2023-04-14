#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 23:35:50 2022

@author: wufuliang
"""

import os
import math
import random
import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB
import time
import networkx as nx
from sklearn.cluster import KMeans # KMeans
from para import *

# =============================================================================
# Number of customers '08', '10', '20', '30', '40', '50', '60', '70', '80'
# =============================================================================
numCus = '08'  # number of customers

# =============================================================================
# Given a tuplelist of edges, find the shortest subtour
# =============================================================================
def subtour(edges):
    unvisited = list(range(ncus))
    cycle = range(ncus+1)  # initial length has 1 more city
    while unvisited:  # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i, j in edges.select(current, '*') if j in unvisited]
        if len(cycle) > len(thiscycle):
            cycle = thiscycle
    return cycle

# =============================================================================
# NetworkX based subtour elimination Cut
# https://networkx.org/documentation/stable/_modules/networkx/algorithms/connectivity/stoerwagner.html
# Stoer, M., & Wagner, F. (1997). A simple min-cut algorithm. Journal of the ACM (JACM), 44(4), 585-591. 
# =============================================================================
def swCut(model,where):
    m._calCou += 1
    if where == GRB.Callback.MIPSOL:
        m._calCouMIPSOL += 1
        Xvals = model.cbGetSolution(model._X)
        usvals = model.cbGetSolution(model._us)
        tour = [ele for ele in model._X.keys() if Xvals[ele] > 0.5]
        
        ### SW minimum cut
        G = nx.Graph()
        G.add_edges_from(tour,weight = 1)
        G.add_edges_from([(0,i) for i in range(1,n+2) if (0,i) not in tour],weight = 0) # Make the connected graph 
        cut_value, partition = nx.stoer_wagner(G)
        if cut_value >= 1:
            for (i,j) in tour:
                # Subgradient cut
                model.cbLazy(model._qs[i,j] >= -1/2*usvals[i,j]**(-3/2)*(model._us[i,j]-model._X[i,j]*usvals[i,j]) + model._X[i,j]*usvals[i,j]**(-1/2))
                # Battery flow
                model.cbLazy(model._ys[i] - model._we[i,j] - mu*model._wb[i,j] - model._wr[i,j] >= model._ys[j] - (1-model._X[i,j]) * Cbat) # (4.60)
                model.cbLazy(model._ys[j] + model._we[i,j] + mu*model._wb[i,j] >= model._ys[i] - (1-model._X[i,j]) * Cbat) # 
                
        else:
            subtour = partition[0]
            subCut = [(i,j) for (i,j) in model._X.keys() if i in subtour and j not in subtour]
            model.cbLazy(gp.quicksum(model._X[i,j] for (i,j) in subCut) >= 1)
            
### Data Read ###
instances = os.listdir('../inst') # instances' file names

n = int(numCus) # Number of customers
ncus = n + 2 # number of customers + 2 ( depot, destination )
    
files = [ele for ele in instances if numCus in ele]

fils, Objs, Objgaps, Caltims, EnerSums, calCous, calMIPCous = [], [], [], [], [], [], [] # objective function and calculation time Dists                                                                                             
for file in files:
    print(files.index(file),file)
    
    if 'HEVTSP_1' in file:
        TB = 3600*1  # Time Budget
    elif 'HEVTSP_2' in file:
        TB = 3600*2  # Time Budget
    else:
        TB = 3600*3  # Time Budget
# =============================================================================
    ### Data Read ###
    corDF = pd.read_csv(open('../inst/' + file), delim_whitespace = True, header = 0) # DataFrame of Coordinates
    
    # Arc list
    arc = []
    for i in range(n+1): # 0 denotes depot
        for j in range(n+1):
            arc.append((i,j))
    for l in range(len(arc)): # 
        if arc[l][1] == 0:
            arc[l] = (arc[l][0],n+1)
     
    # Distance between each nodes
    arcDist, graDist = {}, {} # Arc Dist with Distance, Grade
    for (i,j) in arc:
        if j != n+1:
            arcDist[(i,j)] = round( np.sqrt((corDF['x'][i] - corDF['x'][j])**2 + (corDF['y'][i] - corDF['y'][j])**2), 0 )
        else:
            arcDist[(i,j)] = round( np.sqrt((corDF['x'][i] - corDF['x'][0])**2 + (corDF['y'][i] - corDF['y'][0])**2), 0 )
        
        if (i,j) == (0,n+1) or i == j:
            graDist[i,j] = 0
        else:
            if j != n+1:
                graDist[(i,j)] = math.atan((corDF['z'][j] - corDF['z'][i])/arcDist[i,j])
            else:
                graDist[(i,j)] = math.atan((corDF['z'][0] - corDF['z'][i])/arcDist[i,j])
        
    vUbDist = {} # Upper speed limit
    random.seed(1)
    for (i,j) in arc:
        vUbDist[i,j] = random.uniform(15,19)
        
    EMaxDist = {} # Big M value
    for (i,j) in arc:
        M1 = max( arcDist[i,j]/etaD * (mas*g*np.sin(graDist[i,j]) + Cr*mas*g*np.cos(graDist[i,j])), arcDist[i,j] * etaG * (mas*g*np.sin(graDist[i,j]) + Cr*mas*g*np.cos(graDist[i,j])) ) / 10**6
        M2 = max( arcDist[i,j]/etaD * (mas*g*np.sin(graDist[i,j]) + 0.5*Cd*rho*A*vUbDist[i,j]**2 + Cr*mas*g*np.cos(graDist[i,j])), arcDist[i,j] * etaG * (mas*g*np.sin(graDist[i,j]) +0.5*Cd*rho*A*vUbDist[i,j]**2 + Cr*mas*g*np.cos(graDist[i,j])) ) / 10**6
        EMaxDist[i,j] = max(abs(M1),abs(M2))
    
#     Identify energy consumption
    MlDist,MuDist = {},{}
    for (i,j) in arc:
        MlDist[i,j] = ( mas*g*np.sin(graDist[i,j]) + Cr*mas*g*np.cos(graDist[i,j]) ) * arcDist[i,j]
        MuDist[i,j] = ( mas*g*np.sin(graDist[i,j]) + 0.5*Cd*rho*A*vUbDist[i,j]**2 + Cr*mas*g*np.cos(graDist[i,j]) ) * arcDist[i,j]
    
    ##### Build the model #####
    m = gp.Model()
    
    # Variables for arc level
    arc = [(i,j) for i,j in arc if i!=j]
    X  = m.addVars(arc, vtype=GRB.BINARY, name='x') # x_{ij}
    Xf = m.addVars(arc, vtype=GRB.BINARY, name='xf') # x_{ij}^f
    Xe = m.addVars(arc, vtype=GRB.BINARY, name='xe') # x_{ij}^e
    Xb = m.addVars(arc, vtype=GRB.BINARY, name='xb') # x_{ij}^b
    Xr = m.addVars(arc, vtype=GRB.BINARY, name='xr') # x_{ij}^r
    
    # Fix xr[i,j] for some value
    for (i,j) in arc:
        if MlDist[i,j] > 0:
            m.addConstr(Xr[i,j] == 0)
    
    us  = m.addVars(arc, name='u') # v^2
    qs = m.addVars(arc, lb = 0, name='q') # 
    
    wf = m.addVars(arc, lb = 0, name='wf') # w_{ij}^f 
    we = m.addVars(arc, lb = 0, name='we') # w_{ij}^e 
    wb = m.addVars(arc, lb = 0, name='wb') # w_{ij}^b 
    wr = m.addVars(arc, lb = - GRB.INFINITY, ub = 0,  name='wr')
    Es = m.addVars(arc, lb = -GRB.INFINITY, name='Es')
    ys = m.addVars(np.arange(n+2).tolist(), lb=0, ub=Cbat, name='y') # y_i; (3.14) 
    
    # w Constraints
    m.addConstrs(wf[i,j] <= Xf[i,j] * EMaxDist[i,j] for (i,j) in arc) # (4.2)
    m.addConstrs(wb[i,j] <= Xb[i,j] * EMaxDist[i,j] for (i,j) in arc) # (4.3)
    m.addConstrs(we[i,j] <= Xe[i,j] * EMaxDist[i,j] for (i,j) in arc) # (4.4)
    m.addConstrs(wr[i,j] >= - Xr[i,j] * EMaxDist[i,j] for (i,j) in arc)
    m.addConstrs(wf[i,j] + wb[i,j] + we[i,j] + wr[i,j] >= Es[i,j] for (i,j) in arc)
    
    # Constraint
    m.addConstrs(Es[i,j] >= arcDist[i,j]/etaD * (mas*g*np.sin(graDist[i,j])*X[i,j] + 0.5*Cd*rho*A*us[i,j] + Cr*mas*g*np.cos(graDist[i,j])*X[i,j]) / 10**6  for (i,j) in arc)
    m.addConstrs(Es[i,j] >= arcDist[i,j] * etaG * (mas*g*np.sin(graDist[i,j])*X[i,j] +0.5*Cd*rho*A*us[i,j] + Cr*mas*g*np.cos(graDist[i,j])*X[i,j]) / 10**6 for (i,j) in arc)
    
    # Constraints
    m.addConstrs(X[i,j] == Xf[i,j] + Xe[i,j] + Xb[i,j] + Xr[i,j] for i,j in arc) # (3.8)
    m.addConstrs(gp.quicksum(X[i,j] for (i,j) in arc if i == l and j!=l) == 1 for l in np.arange(n+1)) # (3.9)  
    m.addConstrs(gp.quicksum(X[i,j] for (i,j) in arc if j == l and i!=l) == 1 for l in np.arange(1,n+2)) # (3.10) 
    m.addConstrs( (1 - Xr[i,j]) * EMaxDist[i,j] >= Es[i,j] for (i,j) in arc) 
    m.addConstrs( (Xf[i,j] + Xe[i,j] + Xb[i,j] - 1) * EMaxDist[i,j] <= Es[i,j] for (i,j) in arc) 
    
    # Battery flow
    m.addConstr(ys[0] == Cbat) # (3.16) Right
    
    # Travel time budget
    m.addConstr(gp.quicksum(arcDist[i,j]*qs[i,j] for (i,j) in arc) <= TB ) # (4.45)
    
    # Speed limits
    m.addConstrs(us[i,j] >= vMin**2*X[i,j]  for (i,j) in arc)
    m.addConstrs(us[i,j] <= vUbDist[i,j]**2*X[i,j] for (i,j) in arc)
# =============================================================================
    # Objective function
    obj = m.addVar(lb = -GRB.INFINITY, name='obj') #
    m.addConstr(obj == gp.quicksum(cf * wf[i,j] + ce * we[i,j] + cb * wb[i,j] + ce*wr[i,j] for (i,j) in arc))

# =============================================================================
#     Objective value
    m.setObjective(obj*10**6, GRB.MINIMIZE)
# =============================================================================
    # Branch priority
    km = KMeans(n_clusters = 2)
    arr = np.array([arcDist[ele] for ele in arc if arcDist[ele] > 0]).reshape(-1,1)
    km = km.fit(arr)
    cens = [ele[0] for ele in km.cluster_centers_]
    cenSor = sorted(cens) # sort
    if cenSor[-1]/cenSor[0] < 6:
        for ele in X.keys():
            Xr[ele].setAttr("BranchPriority",100)
            X[ele].setAttr("BranchPriority",80)
            Xf[ele].setAttr("BranchPriority",60)
            Xb[ele].setAttr("BranchPriority",40)
            Xe[ele].setAttr("BranchPriority",20)
    else:
        for ele in X.keys(): 
            label = km.predict(np.array([arcDist[ele]]).reshape(-1,1))[0]
            sor = cenSor.index(cens[label])
            Xr[ele].setAttr("BranchPriority",100+100*sor)
            X[ele].setAttr("BranchPriority",80+100*sor)
            Xf[ele].setAttr("BranchPriority",60+100*sor)
            Xb[ele].setAttr("BranchPriority",40+100*sor)
            Xe[ele].setAttr("BranchPriority",20+100*sor)
            
    m.Params.BranchDir = 1 # Branch up or down
#    m.Params.VarBranch = 1
    m.Params.Threads = 64
# =============================================================================
    m._X,m._qs,m._us=X,qs,us
    m._Xf,m._Xe,m._Xb,m._Xr = Xf,Xe,Xb,Xr
    m._we,m._wb,m._wf,m._wr,m._ys,m._us = we,wb,wf,wr,ys,us
    m._calCou, m._calCouMIPSOL = 0,0
    m.Params.lazyConstraints = 1
    m.setParam("TimeLimit",3600*2) # Limit the solution time; Seconds
    startTime = time.time()
    m.optimize(swCut)
    endTime =  time.time()
    CalTime = endTime - startTime
    print ("--- %s seconds" % (endTime - startTime))
    print('Optimal cost: %g' % m.objVal)
    
    fils.append(file) # File names
    Objs.append(m.objVal)  # Objective values
    Objgaps.append(m.mipgap) # Objective gaps
    Caltims.append(CalTime) # Calculation times
    calCous.append(m._calCou)
    calMIPCous.append(m._calCouMIPSOL)
    
    ResDict = {'Objective':Objs, 'Caltime':Caltims, 'Objgaps':Objgaps, 'Callback': calCous,'CallbackMIP': calMIPCous}
    ResDf = pd.DataFrame(ResDict, index=fils)
    ResDf.to_csv('../result/PHEV' + 'TSPS_'  + numCus + '.csv',index=True,sep=',',encoding='utf_8_sig')
