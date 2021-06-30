# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 12:03:48 2021

@author: User
"""
from Ori_Repres import DCM, DCM2Quaternion, Quaternion_update, DCM2Euler
import numpy as np
p = 0
r = 0
y = 0
dp = 0
dr = 0
dy = 315
V = np.array([[1],[0],[0]])
THE = np.array([[dp],[dr],[dy]])
dt =0.02
[OR, XX] = DCM(p,r,y)
Oq=DCM2Quaternion(OR)

[dR, XX] = DCM(dp,dr,dy)
UR_m1 = dR.dot(OR)
UORI_m1 = DCM2Euler(UR_m1)

[UR_m2, Uq, UORI_m2] = Quaternion_update(THE/dt, Oq, dt)

print(OR)
print(UR_m1)
print(UR_m2)
print(UR_m2-UR_m1)

print(UORI_m1)
print(UORI_m2)
print(UORI_m2-UORI_m1)

print(UR_m2.dot(V))

