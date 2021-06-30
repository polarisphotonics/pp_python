# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:28:22 2021

@author: User
"""
#import tool
import Ref_Elps_Prmt as REP
import Raw_Data_Ld as RDL
import numpy as np
from numpy import sin, pi

#time
hr0 = 12  #<==================Please enter start time (hr-part)
min0 = 13 #<==================Please enter start time (min-part)
sec0 = 46 #<==================Please enter start time (sec-part)
Time = RDL.Raw_data[:,0] + hr0*3600 + min0*60 + sec0 #UTC+8 (unit:sec)

#VBOX date set: Geodetic coordinate/ Ellipsoidal height/ Horizontal & vertical GNSS velocity/ num of satellite
Lat = RDL.Raw_data[:,12] # Latitude 
Lon = RDL.Raw_data[:,13] # Longitude
hei = RDL.Raw_data[:,15] # height
GNSSV = RDL.Raw_data[:,10] # Horizontal GNSS velocity
GNSSV_u = RDL.Raw_data[:,16] # vertical GNSS velocity
Sat_N = RDL.Raw_data[:,14] # num of satellite

#Inertial Sensor

# SRS-200 (single-axis FOG for detecting Wz)
SRS_wz = RDL.Raw_data[:,1] # wz (unit: deg/sec) 

# Sparrow (single-axis FOG produced by Polaris Photonics Ltd for detecting Wz )
PP_wz = RDL.Raw_data[:,2] # original PP_wz (unit: deg/sec)
coeff1 = 0.992463094426231
coeff2 = 5.195494328570396e-04
PP_wz = coeff1*PP_wz + coeff2 # refineed PP_wz (unit: deg/sec)

# Nano-33 (3-axis MEMS gyro for detecting Wx, Wy, Wz)
Nano_wx = RDL.Raw_data[:,3]  # wx (unit: deg/sec)
Nano_wy = RDL.Raw_data[:,4]  # wy (unit: deg/sec)
Nano_wz = RDL.Raw_data[:,5]  # wz (unit: deg/sec)

# Adxl-355 (3-axis MEMS accelerometer for detecting fx, fy, fz)
Adxl_fx  = np.zeros([RDL.num_pen])
Adxl_fy  = np.zeros([RDL.num_pen])
Adxl_fz  = np.zeros([RDL.num_pen]) 

for i in range(RDL.num_pen):
    Lat_temp = RDL.Raw_data[i,12]*(pi/180) # Temporary latitude
    hei_temp = RDL.Raw_data[i,15] # Temporary ellipsoidal height
    g_temp = REP.SF1*(1 + REP.SF2*sin(Lat_temp)**2+ REP.SF3*sin(Lat_temp)**4)\
        + (REP.SF4 + REP.SF5*sin(Lat_temp)**2)*hei_temp + REP.SF6*hei_temp**2 #Temporary normal gravity
    Adxl_fx[i] = RDL.Raw_data[i,7]*g_temp   # fx (unit: m/s^2)
    Adxl_fy[i] = RDL.Raw_data[i,6]*g_temp   # fy (unit: m/s^2)
    Adxl_fz[i] = - RDL.Raw_data[i,8]*g_temp # fz (unit: m/s^2)
