# -*- coding: utf-8 -*-
"""
Created on Mon May 02 10:02:23 2016

@author: ASHWIN KANE
"""
import scipy
import numpy
from scipy.integrate import odeint
import matplotlib.pyplot as plt

ma=1
mb=2
U=100
P=1
L=10
Tain=400
Tbin=300
n=1000
ka=-(U*P/ma)
kb=-(U*P/mb)

def cpa(T):
    f=4000+10*T+0.01*T*T
    return f
def cpadT(T):
    f=4000*T+10*T*T/2+0.01*T*T*T/3
    return f

def cpb(T):
    f=3000+5*T+0.02*T*T
    return f
def cpbdT(T):
    f=3000*T+5*T*T/2+0.02*T*T*T/3
    return f

def dTdx(T,x):
    Ta=T[0]
    Tb=T[1]
    dTadx=(ka/cpa(Ta))*(Ta-Tb)
    Ta=Ta+(dTadx*L/n)
    dTbdx=(kb/cpb(Tb))*(Ta-Tb)
    Tb=Tb+(dTbdx*L/n)
    return scipy.array([dTadx,dTbdx])

Ta0=400
Tb0=200
T0=scipy.array([Ta0,Tb0])
x=scipy.linspace(0,L,n)

T=odeint(dTdx,T0,x)

plt.plot(x,T[:,0],"r-")


































