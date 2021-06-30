# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 21:53:53 2021

@author: User
"""
#import tool
import Init_Orit as IO
from Ori_Repres import ORI2Quaternion, DCM
import numpy as np
from numpy import sqrt

'''' ME's initial value'''
Pitch_PA = IO.Pitch0      # initail value of pitch angle for Plan A
Pitch_PB = IO.Pitch0      # initail value of pitch angle for Plan B
Pitch_PC = IO.Pitch0      # initail value of pitch angle for Plan C

Roll_PA = IO.Roll0       # initail value of roll angle for Plan A
Roll_PB = IO.Roll0       # initail value of roll angle for Plan B
Roll_PC = IO.Roll0       # initail value of roll angle for Plan C

Yaw_PA = 202#158#158 #201 #1.585948249914840e+02 -180#Yaw0_PA;     # initail value of Yaw angle for Plan A 
Yaw_PB = 202#158#158 #1.586566752862757e+02 -180#Yaw0_PB;     # initail value of Yaw angle for Plan B -180
Yaw_PC = 202#158#158 #1.586306747000801e+02 -180#Yaw0_PC;     # initail value of Yaw angle for Plan C -180

heading_PA = Yaw_PA + 180 # initail value of Azimuth angle for Plan A
heading_PB = Yaw_PB + 180 # initail value of Azimuth angle for Plan B
heading_PC = Yaw_PC + 180 # initail value of Azimuth angle for Plan C

''' Orientation (presented by quaternion)'''
[Q_k_PA, DCM_PA] = ORI2Quaternion(Pitch_PA, Roll_PA, Yaw_PA)
[Q_k_PB, DCM_PB] = ORI2Quaternion(Pitch_PB, Roll_PB, Yaw_PB)
[Q_k_PC, DCM_PC] = ORI2Quaternion(Pitch_PC, Roll_PC, Yaw_PC)
ORI_PA = np.array([[Pitch_PA],[Roll_PA],[Yaw_PA]])
ORI_PB = np.array([[Pitch_PB],[Roll_PB],[Yaw_PB]])
ORI_PC = np.array([[Pitch_PC],[Roll_PC],[Yaw_PC]])
[R_PA, R_PAT]= DCM(ORI_PA[0][0], ORI_PA[1][0], ORI_PA[2][0])
[R_PB, R_PBT]= DCM(ORI_PB[0][0], ORI_PB[1][0], ORI_PB[2][0])
[R_PC, R_PCT]= DCM(ORI_PC[0][0], ORI_PC[1][0], ORI_PC[2][0])

''' Velocity'''
VEL_PA = np.array([[0], [0], [0]]) #ve, vn, vu
VEL_PB = np.array([[0], [0], [0]])
VEL_PC = np.array([[0], [0], [0]])

''' Geodetic coordinate'''
POS_PA = np.array([[IO.AVG_Lat], [IO.AVG_Lon], [IO.AVG_hei]])
POS_PB = np.array([[IO.AVG_Lat], [IO.AVG_Lon], [IO.AVG_hei]])
POS_PC = np.array([[IO.AVG_Lat], [IO.AVG_Lon], [IO.AVG_hei]])
'''------------------------------------------------------------------------'''

'''' KF's initial value'''
''' Error State (delta x vector)'''
#items: error of latitude, longitude, height, ve, vn, vu, p, r, y, wx, wy,wz, fx,fy, fz
ERRST_PA = np.zeros([15,1]) #Error State of Plan A
ERRST_PB = np.zeros([15,1]) #Error State of Plan B
ERRST_PC = np.zeros([15,1]) #Error State of Plan C

''' Error State covariance matrix (P matrix)'''
P_PA = np.diag([3,3,9,0.004,0.004,0.004,0.0436,0.0436,0.0436,1.81e-4,8.87e-5,1.75e-5,2.63e-7,2.80e-6,6.99e-7]) #P matrix of Plan A
P_PB = np.diag([3,3,9,0.004,0.004,0.004,0.0436,0.0436,0.0436,1.81e-4,8.87e-5,1.75e-5,2.63e-7,2.80e-6,6.99e-7]) #P matrix of Plan B
P_PC = np.diag([3,3,9,0.004,0.004,0.004,0.0436,0.0436,0.0436,1.81e-4,8.87e-5,1.75e-5,2.63e-7,2.80e-6,6.99e-7]) #P matrix of Plan C

''' Error State noise covariance matrix (Q matrix)'''
Q_PA = np.eye(15,15)
Q_PB = np.eye(15,15)
Q_PC = np.eye(15,15)

''' Gauss-Markov randon process parameters (beta and sigma)'''
btw_PA = np.array([[1/1800], [1/1800], [1/1800]])
btw_PB = np.array([[1/1800], [1/1800], [1/1800]])
btw_PC = np.array([[1/1800], [1/1800], [1/1800]]) 
btf =    np.array([[1/1800], [1/1800], [1/1800]]) 
sgmw_PA = np.array([[1], [1], [1]])
sgmw_PB = np.array([[1], [1], [1]]) 
sgmw_PC = np.array([[1], [1], [1]]) 
sgmf =    np.array([[1], [1], [1]])

''' Jacobian Matrix for Q matrix (G matrix)'''
G_PA = np.diag([IO.STD_Lat**2,IO.STD_Lon**2,IO.STD_hei**2,\
              IO.STD_GNSSV**2,IO.STD_GNSSV**2,IO.STD_GNSSV_u**2,\
              0.01**2,0.01**2,0.01**2,\
              sqrt(2*btw_PA[0][0]*sgmw_PA[0][0]**2),sqrt(2*btw_PA[1][0]*sgmw_PA[1][0]**2),sqrt(2*btw_PA[2][0]*sgmw_PA[2][0]**2),\
              sqrt(2*btf[0][0]*sgmf[0][0]**2),sqrt(2*btf[1][0]*sgmf[1][0]**2),sqrt(2*btf[2][0]*sgmf[2][0]**2)])
G_PB =  np.diag([IO.STD_Lat**2,IO.STD_Lon**2,IO.STD_hei**2,\
              IO.STD_GNSSV**2,IO.STD_GNSSV**2,IO.STD_GNSSV_u**2,\
              0.01**2,0.01**2,0.01**2,\
              sqrt(2*btw_PB[0][0]*sgmw_PB[0][0]**2),sqrt(2*btw_PB[1][0]*sgmw_PB[1][0]**2),sqrt(2*btw_PB[2][0]*sgmw_PB[2][0]**2),\
              sqrt(2*btf[0][0]*sgmf[0][0]**2),sqrt(2*btf[1][0]*sgmf[1][0]**2),sqrt(2*btf[2][0]*sgmf[2][0]**2)])
G_PC = np.diag([IO.STD_Lat**2,IO.STD_Lon**2,IO.STD_hei**2,\
              IO.STD_GNSSV**2,IO.STD_GNSSV**2,IO.STD_GNSSV_u**2,\
              0.01**2,0.01**2,0.01**2,\
              sqrt(2*btw_PC[0][0]*sgmw_PC[0][0]**2),sqrt(2*btw_PC[1][0]*sgmw_PC[1][0]**2),sqrt(2*btw_PC[2][0]*sgmw_PC[2][0]**2),\
              sqrt(2*btf[0][0]*sgmf[0][0]**2),sqrt(2*btf[1][0]*sgmf[1][0]**2),sqrt(2*btf[2][0]*sgmf[2][0]**2)])