# -*- coding: utf-8 -*-
"""
Created on Mon May 02 01:36:50 2016

@author: ASHWIN KANE
"""

import pandas
import scipy

file_excel=pandas.read_excel("ExamProblemData.xlsx","ExamProblemData",header=False,index=False)
t1=file_excel["Col1"]
t1=scipy.array(t1)
print(t1)
c1r=scipy.array(file_excel["Col2"])
print(c1r)
c1p=scipy.array(file_excel["Col3"])
print(c1p)




























































