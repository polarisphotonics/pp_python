# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:15:17 2021

@author: User
"""

from Ori_Repres import DCM, DCM2Quaternion, Quaternion2DCM, DCM2Euler

p = -20
r = +10
y = 200



[R, RT] = DCM(p, r, y)
ORI = DCM2Euler(R)
'''
if R[2][2] >0:
    ORI[1][0] = ORI[1][0]
elif R[2][2] <0 and ORI[1][0] <0:
    ORI[1][0]=ORI[1][0]+180
elif R[2][2] <0 and ORI[1][0] >0:
    ORI[1][0]=ORI[1][0]-180

if R[1][1] >0 and ORI[2][0]>0:
    ORI[2][0] = ORI[2][0]
elif R[1][1] >0 and ORI[2][0]<0:
    ORI[2][0]=ORI[2][0]+360
elif R[1][1] <0:
    ORI[2][0]=ORI[2][0]+180
'''
print(ORI[0][0],ORI[1][0],ORI[2][0])
print(ORI[0][0]-p,ORI[1][0]-r,ORI[2][0]-y)