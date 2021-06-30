# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:47:28 2021

@author: User
"""

import time
tStart = time.time()#計時開始
import numpy as np
from numpy import sqrt, sin, pi, cos
from matplotlib import pyplot as plt
import csv
#**************************** part1: Pre-work****************************
'''Global Variable Declaration'''
import Ref_Elps_Prmt as REP #input geodetic parameter
'''Raw data loading'''
import Raw_Data_Ld as RDL #loading raw data
'''Data Pre-processing'''
import Data_PP as DP #obtain measurement data
'''Initialize Orientation'''
import Init_Orit as IO #calculate the initial value of orientation
'''Initialize parameters for Geodetic Navigation (including ME & KF)'''
import Sph_Nav_init as SNI
'''Import navigation module'''
from State_Est import  Mechanization_Eq
'''Import Kalman Filter module'''
from Kalman_Filter import STM, KF_prediction, KF_correction
'''import Quaternion related module'''
from Ori_Repres import DCM, DCM2Quaternion, Quaternion_update, DCM2Euler, ORI2Quaternion
from AccLevel_Gycomp import Acc_Leveling
'''=====>start point of Geodetic Navigation'''
#Construct zero matrix for saving important data
TJ_PA = np.zeros([RDL.num_pen - 1, 2])   # Trajectory of Plan A
TJ_PB = np.zeros([RDL.num_pen - 1, 2])   # Trajectory of Plan B
TJ_PC = np.zeros([RDL.num_pen - 1, 2])   # Trajectory of Plan C
ORIENTATION_PA = np.zeros([RDL.num_pen - 1, 3])  # Orientation series of Plan A
ORIENTATION_PB = np.zeros([RDL.num_pen - 1, 3])  # Orientation series of Plan B
ORIENTATION_PC = np.zeros([RDL.num_pen - 1, 3])  # Orientation series of Plan C
VELOCITY_PA = np.zeros([RDL.num_pen - 1, 4])  # Velocity series of Plan A
VELOCITY_PB = np.zeros([RDL.num_pen - 1, 4])  # Velocity series of Plan B
VELOCITY_PC = np.zeros([RDL.num_pen - 1, 4])  # Velocity series of Plan C
Length_PA = 0
Length_PB = 0
Length_PC = 0
#initial value for KF
POS_PA = SNI.POS_PA
POS_PB = SNI.POS_PB
POS_PC = SNI.POS_PC
VEL_PA = SNI.VEL_PA
VEL_PB = SNI.VEL_PB
VEL_PC = SNI.VEL_PC
Q_k_PA = SNI.Q_k_PA
Q_k_PB = SNI.Q_k_PB
Q_k_PC = SNI.Q_k_PC
btw_PA = SNI.btw_PA
btw_PB = SNI.btw_PB
btw_PC = SNI.btw_PC
btf = SNI.btf
ERRST_PA = SNI.ERRST_PA
ERRST_PB = SNI.ERRST_PB
ERRST_PC = SNI.ERRST_PC
G_PA = SNI.G_PA
G_PB = SNI.G_PB
G_PC = SNI.G_PC
Q_PA = SNI.Q_PA
Q_PB = SNI.Q_PB
Q_PC = SNI.Q_PC
P_PA = SNI.P_PA
P_PB = SNI.P_PB
P_PC = SNI.P_PC
Length_PA = 0
Length_PB = 0
Length_PC = 0
AZ_PA = SNI.heading_PA
AZ_PB = SNI.heading_PB
AZ_PC = SNI.heading_PC
ORI_PA = SNI.ORI_PA
ORI_PB = SNI.ORI_PB
ORI_PC = SNI.ORI_PC
R_PA = SNI.R_PA
R_PB = SNI.R_PB
R_PC = SNI.R_PC
ken =[]
ken1=[]
ken2=[]
for i in range(RDL.num_pen - 1):
    dt = DP.Time[i+1] - DP.Time[i]
    F_k = np.array([[DP.Adxl_fx[i]], [DP.Adxl_fy[i]], [DP.Adxl_fz[i]]]) #-(ERRST_PA[12:15]+ERRST_PB[12:15]+ERRST_PC[12:15])/3
    W_k_PA = np.array([[DP.Nano_wx[i]], [DP.Nano_wy[i]], [DP.SRS_wz[i]]]) #-ERRST_PA[9:12]
    W_k_PB = np.array([[DP.Nano_wx[i]], [DP.Nano_wy[i]], [DP.PP_wz[i]]])  #-ERRST_PB[9:12]
    W_k_PC = np.array([[DP.Nano_wx[i]], [DP.Nano_wy[i]], [DP.Nano_wz[i]]]) #-ERRST_PC[9:12]
    [Pitch_w, Roll_w] = Acc_Leveling(F_k[0][0], F_k[1][0], F_k[2][0], 9.8)
    ORI_k = np.array([[Pitch_w],[Roll_w]])
    ken = np.append(ken,Pitch_w)
    ken1 = np.append(ken1,Roll_w)
    #****************************Part2: State Estimation by ME****************************
    '''POS_k1, VEL_k1, ORI, Q_k1, R_k1, omg_LBB = Mechanization_Eq(POS_k, VEL_k, Q_k, W_k, F_k, dt)'''
    
    [POS_PA, VEL_PA, R_PA, ORI_PA, omg_LBB_PA] = Mechanization_Eq(POS_PA, VEL_PA, R_PA, W_k_PA, F_k, dt)
    [POS_PB, VEL_PB, R_PB, ORI_PB, omg_LBB_PB] = Mechanization_Eq(POS_PB, VEL_PB, R_PB, W_k_PB, F_k, dt)
    [POS_PC, VEL_PC, R_PC, ORI_PC, omg_LBB_PC] = Mechanization_Eq(POS_PC, VEL_PC, R_PC, W_k_PC, F_k, dt)
    
    #d its covariance for KF'''
    # input external measurement Zk (from GNSS) & Enter New Measurement Covariance
    
    if DP.Sat_N[i] == 0 and i >70000 and i <75000: #因應隧道速限80kmh，如果無衛星，速度直接當定速
       DP.GNSSV[i]=16.6666666667
    
    #GNSS Velocity of Plan A, B, C
    VELg_PA = np.array([[DP.GNSSV[i]*sin((ORI_PA[2][0])*pi/180)], [DP.GNSSV[i]*cos((ORI_PA[2][0])*pi/180)], [DP.GNSSV_u[i]]])
    VELg_PB = np.array([[DP.GNSSV[i]*sin((ORI_PB[2][0])*pi/180)], [DP.GNSSV[i]*cos((ORI_PB[2][0])*pi/180)], [DP.GNSSV_u[i]]])
    VELg_PC = np.array([[DP.GNSSV[i]*sin((ORI_PC[2][0])*pi/180)], [DP.GNSSV[i]*cos((ORI_PC[2][0])*pi/180)], [DP.GNSSV_u[i]]])
    
    #GNSS position
    POSg = np.array([[DP.Lat[i]], [DP.Lon[i]], [DP.hei[i]]])
    
    #Zk for KF correction
    Zk_PA = -np.vstack([POSg,VELg_PA,ORI_k]) + np.vstack([POS_PA,VEL_PA,ORI_PA[0:2]])
    Zk_PB = -np.vstack([POSg,VELg_PB,ORI_k]) + np.vstack([POS_PB,VEL_PB,ORI_PB[0:2]])
    Zk_PC = -np.vstack([POSg,VELg_PC,ORI_k]) + np.vstack([POS_PC,VEL_PC,ORI_PC[0:2]])
    
    #Rk for KF correction
    
   
    '''
    if i>70000:
        if DP.Sat_N[i] < 3:
            Rk = np.diag([3000000,3000000,9,0.0004,0.0004,0.0004,0.01,0.01])*0.0001
        elif DP.Sat_N[i] >= 3 and DP.Sat_N[i] < 5:
            Rk = np.diag([0.1*3000000,0.1*3000000,9,0.0004,0.0004,0.0004,0.01,0.01])*0.0001
        elif DP.Sat_N[i] >= 5:
            Rk = np.diag([0.1*3000000,0.1*3000000,9,0.0004,0.0004,0.0004,0.01,0.01])*0.00001
    elif i<70000:
        if DP.Sat_N[i] < 3:
            Rk = np.diag([3,3,9,0.0004,0.0004,0.0004,0.01,0.01])*0.0001
        elif DP.Sat_N[i] >= 3 and DP.Sat_N[i] < 5:
            Rk = np.diag([0.1*3,0.1*3,9,0.0004,0.0004,0.0004,0.01,0.01])*0.0001
        elif DP.Sat_N[i] >= 5:
            Rk = np.diag([0.1*3,0.1*3,9,0.0004,0.0004,0.0004,0.01,0.01])*0.00001
    '''     
    Rk = np.diag([3,3,9,0.0004,0.0004,0.0004,0.01,0.01])*0.001
    
    UPD_freq = 3#3 # KF correction period
    SAMP_freq = 85#85 # sampleing frequency (85HZ)
    
    
        
    '''Kalman Filter processing'''
    if (i==0 or i%SAMP_freq*UPD_freq==0):# The period of KF correction is UPD_freq(sec)
        LO_PA = STM(POS_PA, VEL_PA, R_PA, F_k, btw_PA, btf, dt)
        LO_PB = STM(POS_PB, VEL_PB, R_PB, F_k, btw_PB, btf, dt)
        LO_PC = STM(POS_PC, VEL_PC, R_PC, F_k, btw_PC, btf, dt)
    
    # KF prediction
        [ERRST_PA, P_PA] = KF_prediction (ERRST_PA, LO_PA, G_PA, Q_PA, P_PA, dt)
        [ERRST_PB, P_PB] = KF_prediction (ERRST_PB, LO_PB, G_PB, Q_PB, P_PB, dt)
        [ERRST_PC, P_PC] = KF_prediction (ERRST_PC, LO_PC, G_PC, Q_PC, P_PC, dt)    
    
    #****************************Part3: Optimal Error State Estimation****************************
        [ERRST_PA, P_PA] = KF_correction (ERRST_PA, P_PA, Zk_PA, Rk)
        [ERRST_PB, P_PB] = KF_correction (ERRST_PB, P_PB, Zk_PB, Rk)
        [ERRST_PC, P_PC] = KF_correction (ERRST_PC, P_PC, Zk_PC, Rk)
        #print(ERRST_PA)    
        # position compensation
        POS_PA = POS_PA - ERRST_PA[0:3]
        POS_PB = POS_PB - ERRST_PB[0:3]
        POS_PC = POS_PC - ERRST_PC[0:3]
        
        # velocity compensation
        VEL_PA = VEL_PA - ERRST_PA[3:6]
        VEL_PB = VEL_PB - ERRST_PB[3:6]
        VEL_PC = VEL_PC - ERRST_PC[3:6]
        #****************************Part4: Optimal State Estimation****************************
        
        print(ORI_PA)
        print(ERRST_PA)
        ORI_PA[0:2] = ORI_PA[0:2] - np.array([[ERRST_PA[6][0]],[ERRST_PA[7][0]]])
        [R_PA,XXX] = DCM(ORI_PA[0][0], ORI_PA[1][0], ORI_PA[2][0])
        
        ORI_PB[0:2] = ORI_PB[0:2] - np.array([[ERRST_PB[6][0]],[ERRST_PB[7][0]]])
        [R_PB,XXX] = DCM(ORI_PB[0][0], ORI_PB[1][0], ORI_PB[2][0])
        
        ORI_PC[0:2] = ORI_PC[0:2] - np.array([[ERRST_PC[6][0]],[ERRST_PC[7][0]]])
        [R_PC,XXX] = DCM(ORI_PC[0][0], ORI_PC[1][0], ORI_PC[2][0])
        '''
        ORI_PA = ORI_PA + np.array([[ERRST_PA[6][0]],[ERRST_PA[7][0]],[ERRST_PA[8][0]]])
        [R_PA,XXX] = DCM(ORI_PA[0][0], ORI_PA[1][0], ORI_PA[2][0])
        '''
        '''
        ORI_PA = ORI_PA -np.array([[ERRST_PA[6][0]],[ERRST_PA[7][0]],[ERRST_PA[8][0]]])
        ORI_PB = ORI_PB -np.array([[ERRST_PB[6][0]],[ERRST_PB[7][0]],[ERRST_PB[8][0]]])
        ORI_PC = ORI_PC -np.array([[ERRST_PC[6][0]],[ERRST_PC[7][0]],[ERRST_PC[8][0]]])
        [R_PA,XXX] = DCM(ORI_PA[0][0], ORI_PA[1][0], ORI_PA[2][0])
        [R_PB,XXX] = DCM(ORI_PB[0][0], ORI_PB[1][0], ORI_PB[2][0])
        [R_PC,XXX] = DCM(ORI_PC[0][0], ORI_PC[1][0], ORI_PC[2][0])
        '''
    
    #orientation compensation
    #print(ORI_PA)
    
    #print(ERRST_PA[6:9])
    '''
    #print(ORI_PA)
    
    #print(ORI_PA)
    S_PA = np.array([[0, -ERRST_PA[8][0], ERRST_PA[7][0]],
                     [ERRST_PA[8][0], 0, -ERRST_PA[6][0]],
                     [-ERRST_PA[7][0], ERRST_PA[6][0], 0]])
    R_PA = R_PA.dot(np.eye(3,3)+S_PA*pi/180)
    
    '''
    #****************************Part5: Export important information****************************
    
    # path length
    Length_PA = Length_PA + sqrt((VEL_PA[0][0])**2+(VEL_PA[1][0])**2)*dt
    Length_PB = Length_PB + sqrt((VEL_PB[0][0])**2+(VEL_PB[1][0])**2)*dt
    Length_PC = Length_PC + sqrt((VEL_PC[0][0])**2+(VEL_PC[1][0])**2)*dt
    
    # Trajectory
    TJ_PA[i,:] = np.array([POS_PA[1][0], POS_PA[0][0]])
    TJ_PB[i,:] = np.array([POS_PB[1][0], POS_PB[0][0]])
    TJ_PC[i,:] = np.array([POS_PC[1][0], POS_PC[0][0]])
    
    # Velocity
    VELOCITY_PA[i,:] = np.array([VEL_PA[0][0], VEL_PA[1][0], VEL_PA[2][0],\
                                 sqrt((VEL_PA[0][0])**2+(VEL_PA[1][0])**2)])
    VELOCITY_PB[i,:] = np.array([VEL_PB[0][0], VEL_PB[1][0], VEL_PB[2][0],\
                                 sqrt((VEL_PB[0][0])**2+(VEL_PB[1][0])**2)])
    VELOCITY_PC[i,:] = np.array([VEL_PC[0][0], VEL_PC[1][0], VEL_PC[2][0],\
                                 sqrt((VEL_PC[0][0])**2+(VEL_PC[1][0])**2)])
    
    # Orientation
    ORIENTATION_PA[i,:] = np.array([ORI_PA[0][0], ORI_PA[1][0], ORI_PA[2][0]])
    ORIENTATION_PB[i,:] = np.array([ORI_PB[0][0], ORI_PB[1][0], ORI_PB[2][0]])
    ORIENTATION_PC[i,:] = np.array([ORI_PC[0][0], ORI_PC[1][0], ORI_PC[2][0]])
#****************************Part6: Cartographic analysis****************************
'''figure1 for plot TJ'''
plt.rcParams["figure.figsize"] = (10, 10)
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(TJ_PA[:,0], TJ_PA[:,1],'r-'\
        ,TJ_PB[:,0], TJ_PB[:,1],'b-'\
        ,TJ_PC[:,0], TJ_PC[:,1],'k-'\
        ,DP.Lon, DP.Lat,'g-')
plt.legend(labels=['TJ_PA','TJ_PB','TJ_PC','VBOX'],loc='best')
ax.set_aspect('equal')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Trajectory of the three plans')
plt.grid()

'''figure2 for plot VEL'''
plt.rcParams["figure.figsize"] = (20, 20)
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ironman_grid = plt.GridSpec(5, 3, wspace=0.2, hspace=0.5)
plt.subplot(ironman_grid[0,0])
plt.plot(DP.Time[0:RDL.num_pen-1],VELOCITY_PA[:,0]*3.6,'r-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PB[:,0]*3.6,'b-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PC[:,0]*3.6,'k-')
plt.legend(labels=['PA','PB','PC'],loc='best') 
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Velocity (km/h)')
plt.title('Ve')
plt.grid()

plt.subplot(ironman_grid[0,1])
plt.plot(DP.Time[0:RDL.num_pen-1],VELOCITY_PA[:,1]*3.6,'r-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PB[:,1]*3.6,'b-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PC[:,1]*3.6,'k-')
plt.legend(labels=['PA','PB','PC'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Velocity (km/h)') 
plt.title('Vn')
plt.grid()

plt.subplot(ironman_grid[0,2])
plt.plot(DP.Time[0:RDL.num_pen-1],VELOCITY_PA[:,2]*3.6,'r-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PB[:,2]*3.6,'b-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PC[:,2]*3.6,'k-')
plt.legend(labels=['PA','PB','PC'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Velocity (km/h)') 
plt.title('Vu')
plt.grid()

plt.subplot(ironman_grid[1,:])
plt.plot(DP.Time[0:RDL.num_pen-1],VELOCITY_PA[:,3]*3.6,'r-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PB[:,3]*3.6,'b-',\
         DP.Time[0:RDL.num_pen-1],VELOCITY_PC[:,3]*3.6,'k-',\
         DP.Time[0:RDL.num_pen-1],DP.GNSSV[0:RDL.num_pen-1]*3.6,'g-')
plt.legend(labels=['PA','PB','PC','GNSSV'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Velocity (km/h)') 
plt.title('Horizontal V')
plt.grid()

plt.subplot(ironman_grid[2,:])
plt.scatter(DP.Time[0:RDL.num_pen-1],DP.Sat_N[0:RDL.num_pen-1])
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Sat_N') 
plt.title('number of Sat')
plt.grid()

plt.subplot(ironman_grid[3,:])
plt.plot(DP.Time[0:RDL.num_pen-1],TJ_PA[:,0],'r-')
plt.plot(DP.Time[0:RDL.num_pen-1],TJ_PB[:,0],'b-')
plt.plot(DP.Time[0:RDL.num_pen-1],TJ_PC[:,0],'k-')
plt.plot(DP.Time[0:RDL.num_pen-1],DP.Lon[0:RDL.num_pen-1],'g-')
plt.legend(labels=['TJ_PA','TJ_PB','TJ_PC','VBOX'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Lon (deg)') 
plt.title('Longitude')
plt.grid()

plt.subplot(ironman_grid[4,:])
plt.plot(DP.Time[0:RDL.num_pen-1],TJ_PA[:,1],'r-')
plt.plot(DP.Time[0:RDL.num_pen-1],TJ_PB[:,1],'b-')
plt.plot(DP.Time[0:RDL.num_pen-1],TJ_PC[:,1],'k-')
plt.plot(DP.Time[0:RDL.num_pen-1],DP.Lat[0:RDL.num_pen-1],'g-')
plt.legend(labels=['TJ_PA','TJ_PB','TJ_PC','VBOX'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('Lat (deg)') 
plt.title('Latitude')
plt.grid()

'''figure3 for plot ORI'''
plt.rcParams["figure.figsize"] = (20, 10)
fig3 = plt.figure()
ax = fig3.add_subplot(111)
ironman_grid = plt.GridSpec(1, 3, wspace=0.2, hspace=0.3)
plt.subplot(ironman_grid[0,0])
plt.plot(DP.Time[0:RDL.num_pen-1],ORIENTATION_PA[:,0],'r-',\
         DP.Time[0:RDL.num_pen-1],ORIENTATION_PB[:,0],'b-',\
         DP.Time[0:RDL.num_pen-1],ORIENTATION_PC[:,0],'k-')
plt.legend(labels=['PA','PB','PC'],loc='best') 
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('ORIENTATION (ded)')
plt.title('Pitch')
plt.grid()

plt.subplot(ironman_grid[0,1])
plt.plot(DP.Time[0:RDL.num_pen-1],ORIENTATION_PA[:,1],'r-',\
         DP.Time[0:RDL.num_pen-1],ORIENTATION_PB[:,1],'b-',\
         DP.Time[0:RDL.num_pen-1],ORIENTATION_PC[:,1],'k-')
plt.legend(labels=['PA','PB','PC'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('ORIENTATION (deg)') 
plt.title('Roll')
plt.grid()

plt.subplot(ironman_grid[0,2])
plt.plot(DP.Time[0:RDL.num_pen-1],ORIENTATION_PA[:,2],'r-',\
         DP.Time[0:RDL.num_pen-1],ORIENTATION_PB[:,2],'b-',\
         DP.Time[0:RDL.num_pen-1],ORIENTATION_PC[:,2],'k-')
plt.legend(labels=['PA','PB','PC'],loc='best')
plt.xlabel('Time (UTC+8, sec)')
plt.ylabel('ORIENTATION (deg)') 
plt.title('Yaw')
plt.grid()

#****************************Part7: Text information****************************
tEnd = time.time()#計時結束
print( "It cost %f sec \n" % (tEnd - tStart))#會自動做進位
print("Path length of TJ: PA = %7.4f (km)\n" % (Length_PA/1000))
print("Path length of TJ: PB = %7.4f (km)\n" % (Length_PB/1000))
print("Path length of TJ: PC = %7.4f (km)\n" % (Length_PC/1000))

with open('TJ'+'.csv', 'w', newline='') as csvfile:
  # 建立 CSV 檔寫入器
  writer = csv.writer(csvfile)

  # 寫入一列資料
  writer.writerow(['TJ_PA_lon', 'TJ_PA_lat','TJ_PB_lon', 'TJ_PB_lat','TJ_PC_lon', 'TJ_PC_lat'])

  # 寫入另外幾列資料
  for j in range(RDL.num_pen - 1):
      writer.writerow([TJ_PA[j,0], TJ_PA[j,1], TJ_PB[j,0], TJ_PB[j,1],TJ_PC[j,0], TJ_PC[j,1]])

with open('ORI'+'.csv', 'w', newline='') as csvfile:
  # 建立 CSV 檔寫入器
  writer = csv.writer(csvfile)

  # 寫入一列資料
  writer.writerow(['ORI_PA_p', 'ORI_PA_r','ORI_PA_A', 'ORI_PB_p', 'ORI_PB_r','ORI_PB_A','ORI_PC_p', 'ORI_PC_r','ORI_PC_A'])

  # 寫入另外幾列資料
  for j in range(RDL.num_pen - 1):
      writer.writerow([ORIENTATION_PA[j,0], ORIENTATION_PA[j,1], ORIENTATION_PA[j,2],\
                       ORIENTATION_PB[j,0], ORIENTATION_PB[j,1], ORIENTATION_PB[j,2],\
                       ORIENTATION_PC[j,0], ORIENTATION_PC[j,1], ORIENTATION_PC[j,2]])
      