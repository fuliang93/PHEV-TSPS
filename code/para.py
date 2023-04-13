#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 14:37:24 2023

@author: wufuliang
"""

tMax = 10**3 # Big constant
etaD = 0.45 # driventrain Efficiency
etaG = 0.15 # regeneration efficiency
mas = 6350 # kg
g = 9.81 # m/s^2
Cd = 0.7 # Coefficient of aerodynamic drag
rho = 1.2041 # kg/m^3
A = 3.912 # m^2
Cr = 0.01 # Coefficient of rolling resistance
Pacc = 0 # power of Accessory
Cbat = 14.4 * 3600000 / 10**6  # energy limit of Battery 14.4 kw*h J
mu = 0.5 # battery split of energy
cf, cb, ce = 1, 0.7, 0.5

vMin,vMax = 3,19 # Minimum, maximum speeds

epsilon = 60000/ 10**6  # (j/s) charging rate 60 kwh/h 


