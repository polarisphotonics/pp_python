# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 09:58:54 2021

@author: User
"""
import Ref_Elps_Prmt as REP
from numpy import sin, pi, arctan, cos, sqrt
from Ori_Repres import DCM, DCM2Euler, DCM2Quaternion, Quaternion_update, Quaternion2DCM
import numpy as np

def Mechanization_Eq(POS_k, VEL_k, R_k, W_k, F_k, dt):
    '''input data'''
    lat_k = POS_k[0][0]
    lon_k = POS_k[1][0]
    hei_k = POS_k[2][0]
    ve_k = VEL_k[0][0]
    vn_k = VEL_k[1][0]
    vu_k = VEL_k[2][0]
    wx_k = W_k[0][0]
    wy_k = W_k[1][0]
    wz_k = W_k[2][0]
    fx_k = F_k[0][0]
    fy_k = F_k[1][0]
    fz_k = F_k[2][0]
    
    '''Attitude Computation'''
    omg_IBB = np.array([[wx_k], [wy_k], [wz_k]]) # observation of rotation rate
    RN = (REP.Major_R)/(1-(REP.Eccen*sin(lat_k*pi/180))**2)**0.5 # normal radius
    RM = (REP.Major_R*(1-REP.Eccen**2))/(1-(REP.Eccen*sin(lat_k*pi/180))**2)**1.5 # meridian radius
    omg_ELL = (180/pi)*np.array([[-vn_k/(RM+hei_k)],
                                 [ve_k/(RN+hei_k)],
                                 [ve_k*arctan(lat_k*pi/180)/((RN+hei_k))]]) # body motion
    #R_k = Quaternion2DCM(Q_k)
    #[R_k,R_kT] = DCM(ORI_k[0][0],ORI_k[1][0] , ORI_k[2][0])
    omg_LBB = omg_IBB - (R_k.T).dot(omg_ELL)#####################
    #print(R_k)
    #Q_k =DCM2Quaternion(R_k)
    #[R_k1, Q_k1, ORI] = Quaternion_update(omg_LBB, Q_k, dt)
    #print(R_k1)
    #[THE, THET] = DCM(omg_LBB[0][0]*dt,omg_LBB[1][0]*dt,omg_LBB[2][0]*dt)
    
    THE = np.array([[0, -omg_LBB[2][0]*dt, omg_LBB[1][0]*dt],
                    [omg_LBB[2][0]*dt, 0, -omg_LBB[0][0]*dt],
                    [-omg_LBB[1][0]*dt, omg_LBB[0][0]*dt, 0]])
    R_k1 = R_k.dot(np.eye(3,3)+THE*pi/180)
    ORI = DCM2Euler(R_k1)
    #print(ORI)
    '''coordinate transformation for specific force'''
    f_B = np.array([[fx_k], [fy_k], [fz_k]])
    OMG_ELL = np.array([[0, -omg_ELL[2][0],  omg_ELL[1][0]],
                        [omg_ELL[2][0],  0, -omg_ELL[0][0]],
                        [-omg_ELL[1][0], omg_ELL[0][0], 0]])
    VEL_L = np.array([[ve_k], [vn_k], [vu_k]])
    f_L = R_k1.dot(f_B) - OMG_ELL.dot(VEL_L)
    ve_k1 = ve_k + f_L[0][0]*dt
    vn_k1 = vn_k + f_L[1][0]*dt
    vu_k1 = vu_k + f_L[2][0]*dt #####################原本是負號
    VEL_k1 = np.array([[ve_k1], [vn_k1], [vu_k1]])
    
    hei_k1 = hei_k + 0.5*(vu_k1+vu_k)*dt
    lat_k1 = lat_k + (180/pi)*((0.5*(vn_k1+vn_k)*dt)/(RN+hei_k1))
    lon_k1 = lon_k + (180/pi)*((0.5*(ve_k1+ve_k)*dt)/((RN+hei_k1)*cos(lat_k1*pi/180)))
    POS_k1 = np.array([[lat_k1], [lon_k1], [hei_k1]])
    return POS_k1, VEL_k1, R_k1, ORI, omg_LBB