# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 16:45:32 2021

@author: User
"""
#import tool
import Ref_Elps_Prmt as REP
import Raw_Data_Ld as RDL
import Data_PP as DP
import numpy as np
from numpy import pi, sin
from AccLevel_Gycomp import Acc_Leveling, Gyrocompassing

''' Moving Tetection'''
for i  in range(RDL.num_pen):
    if abs(DP.GNSSV[i]) > 0.01:
        Stamp_Start = i
        break
print('The start time of vehicle = %7.4f sec\n' % DP.Time[Stamp_Start]) #ouput start time

'''Data Truncation'''
back = 10
# SRS-200 (single-axis FOG for detecting Wz)
AVG_SRS_wz = np.mean(DP.SRS_wz[0:Stamp_Start-back])   # Average of Wz before driving
STD_SRS_wz =  np.std(DP.SRS_wz[0:Stamp_Start-back])   # Standard deviation of Wz before driving

# Sparrow (single-axis FOG produced by Polaris Photonics Ltd for detecting Wz )
AVG_PP_wz = np.mean(DP.PP_wz[0:Stamp_Start-back])     # Average of Wz before driving
STD_PP_wz =  np.std(DP.PP_wz[0:Stamp_Start-back])     # Standard deviation of Wz before driving

# Nano-33 (3-axis MEMS gyro for detecting Wx, Wy, Wz)
AVG_Nano_wx = np.mean(DP.Nano_wx[0:Stamp_Start-back]) # Average of Wx before driving
AVG_Nano_wy = np.mean(DP.Nano_wy[0:Stamp_Start-back]) # Average of Wy before driving
AVG_Nano_wz = np.mean(DP.Nano_wz[0:Stamp_Start-back]) # Average of Wz before driving
STD_Nano_wx =  np.std(DP.Nano_wx[0:Stamp_Start-back]) # Standard deviation of Wx before driving
STD_Nano_wy =  np.std(DP.Nano_wy[0:Stamp_Start-back]) # Standard deviation of Wy before driving
STD_Nano_wz =  np.std(DP.Nano_wz[0:Stamp_Start-back]) # Standard deviation of Wz before driving

# Adxl-355 (3-axis MEMS accelerometer for detecting fx, fy, fz)
AVG_fx = np.mean(DP.Adxl_fx[0:Stamp_Start-back]) # Average of fx before driving
AVG_fy = np.mean(DP.Adxl_fy[0:Stamp_Start-back]) # Average of fy before driving
AVG_fz = np.mean(DP.Adxl_fz[0:Stamp_Start-back]) # Average of fz before driving
STD_fx =  np.std(DP.Adxl_fx[0:Stamp_Start-back]) # Standard deviation of fx before driving
STD_fy =  np.std(DP.Adxl_fy[0:Stamp_Start-back]) # Standard deviation of fy before driving
STD_fz =  np.std(DP.Adxl_fz[0:Stamp_Start-back]) # Standard deviation of fz before driving

# GNSS information
AVG_Lat = np.mean(DP.Lat[0:Stamp_Start-back])        # Average of Latitude before driving
STD_Lat = np.std(DP.Lat[0:Stamp_Start-back])         # Standard deviation of Latitude before driving
AVG_Lon = np.mean(DP.Lon[0:Stamp_Start-back])        # Average of Longitude before driving
STD_Lon = np.std(DP.Lon[0:Stamp_Start-back])         # Standard deviation of Longitude before driving
AVG_hei = np.mean(DP.hei[0:Stamp_Start-back])        # Average of height before driving
STD_hei = np.std(DP.hei[0:Stamp_Start-back])          # Standard deviation of height before driving
AVG_GNSSV = np.mean(DP.GNSSV[0:Stamp_Start-back])    # Average of horizontal velocity before driving
STD_GNSSV = np.std(DP.GNSSV[0:Stamp_Start-back])     # Standard deviation of horizontal velocity before driving
AVG_GNSSV_u = np.mean(DP.GNSSV_u[0:Stamp_Start-back])# Average of vertical velocity before driving
STD_GNSSV_u = np.std(DP.GNSSV_u[0:Stamp_Start-back]) # Standard deviation of vertical velocity before driving

#Illustraction for Plan A, B, C
'''
Only goal: Trajectory Computation
It has the three data combinations:
Plan A (abbr.PA): Adxl-355 (fx, fy, fz) + Nano-33 (Wx, Wy) + SRS-200 (Wz)
Plan B (abbr.PB): Adxl-355 (fx, fy, fz) + Nano-33 (Wx, Wy) + Sparrow (Wz)
Plan C (abbr.PC): Adxl-355 (fx, fy, fz) + Nano-33 (Wx, Wy, Wz)
'''
PHI0 = DP.Lat[0]*pi/180 # Latitude of origin
HEI0 = DP.hei[0] # Ellisoidal height of origin
g_start = REP.SF1*(1 + REP.SF2*sin(PHI0*pi/180)**2+ REP.SF3*sin(PHI0*pi/180)**4) + \
          (REP.SF4 + REP.SF5*sin(PHI0*pi/180)**2)*HEI0 + REP.SF6*(HEI0)**2 # normal gravity of origin
[Pitch0, Roll0] = Acc_Leveling(AVG_fx, AVG_fy, AVG_fz, g_start)

Yaw0_PA = Gyrocompassing(AVG_Nano_wx, AVG_Nano_wy, AVG_SRS_wz, Pitch0, Roll0, PHI0)
Yaw0_PB = Gyrocompassing(AVG_Nano_wx, AVG_Nano_wy, AVG_PP_wz, Pitch0, Roll0, PHI0)
Yaw0_PC = Gyrocompassing(AVG_Nano_wx, AVG_Nano_wy, AVG_Nano_wz, Pitch0, Roll0, PHI0)