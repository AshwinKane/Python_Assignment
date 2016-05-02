# -*- coding: utf-8 -*-
"""
Created on Mon May 02 10:49:34 2016

@author: ASHWIN KANE
"""
n=10
def f(x):
    f=1/(x**n)
    return f
    
def integrate(x0,xn,nr):
    w=(float(xn)-float(x0))/nr
    s=0
    for i in range(nr):
        h=f(x0+i*w)
        area=h*w
        s+=area
    return s

y=integrate(50,2,100000)
print(y)



















































