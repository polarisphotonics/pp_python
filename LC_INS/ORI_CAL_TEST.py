# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 10:15:17 2021

@author: User
"""

from Ori_Repres import DCM, DCM2Quaternion, Quaternion2DCM, DCM2Euler

p = -20
r = +10
y = 315
'''
if y>0 and y<90:
    yp=y
elif y>90 and y<180:
    yp=y-180
elif y>180 and y<270:
    yp=y-180
elif y>270 and y<360:
    yp=y-360
'''

[R, RT] = DCM(p, r, y)
print(R)
q = DCM2Quaternion(R)
print(q)
Rp = Quaternion2DCM(q)
print(Rp)
ORI = DCM2Euler(Rp)

print(ORI[0][0],ORI[1][0],ORI[2][0])
print(ORI[0][0]-p,ORI[1][0]-r,ORI[2][0]-y)