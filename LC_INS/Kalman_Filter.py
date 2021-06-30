# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:34:58 2021

@author: User
"""
import Ref_Elps_Prmt as REP
from numpy import sqrt, sin, pi, cos, tan
import numpy as np
from Ori_Repres import DCM

def STM(POS, VEL, R, ACC, betaw, betaf, dt):
    phi = POS[0][0]
    hei = POS[2][0]
    #ve = VEL[0][0]
    #vn = VEL[1][0]
    #vu = VEL[2][0]
    fx = ACC[0][0]
    fy = ACC[1][0]
    fz = ACC[2][0]
    #omg_iee = 0 # for 0422路試資料使用，因為移除柯氏力訊號，所以跟柯氏力相關的term皆為零，可以透過給omg_iee = 0達成此目的(有確認過公式可以這樣改)
    RN = (REP.Major_R)/(1-(REP.Eccen*sin(phi*pi/180))**2)**0.5 # normal radius
    RM = (REP.Major_R*(1-REP.Eccen**2))/(1-(REP.Eccen*sin(phi*pi/180))**2)**1.5 # meridian radius
    ZM = np.zeros([3,3])
    fe = R[0][0]*fx + R[0][1]*fy +R[0][2]*fz
    fn = R[1][0]*fx + R[1][1]*fy +R[1][2]*fz
    fu = R[2][0]*fx + R[2][1]*fy +R[2][2]*fz
    '''F11~F15'''
    F11 = ZM
    F12 = np.array([[0, 1/(RM+hei), 0], [1/((RN+hei)*(cos(phi*pi/180))), 0, 0], [0, 0, 1]])
    F13 = ZM
    F14 = ZM
    F15 = ZM
    F1 = np.hstack([F11, F12, F13, F14, F15])
    '''F21~F25'''
    F21 = ZM
    F22 = ZM
    F23 = np.array([[0, fu, -fn], [-fu, 0, fe], [fn, -fe, 0]])
    F24 = ZM
    F25 = R
    F2 = np.hstack([F21, F22, F23, F24, F25])
    '''F31~F35'''
    F31 = ZM
    F32 = np.array([[0,  1/(RM+hei), 0],
                    [-1/(RN+hei), 0, 0],
                    [-tan(phi*pi/180)/(RN+hei), 0, 0]])
    F33 = ZM
    F34 = R
    F35 = ZM
    F3 = np.hstack([F31, F32, F33, F34, F35])
    '''F41~F45'''
    F41 = ZM
    F42 = ZM
    F43 = ZM
    F44 = np.diag([-betaw[0][0],-betaw[1][0],-betaw[2][0]]) #陀螺儀隨機模式矩陣
    F45 = ZM
    F4 = np.hstack([F41, F42, F43, F44, F45])
    '''F51~F55'''
    F51 = ZM
    F52 = ZM
    F53 = ZM
    F54 = ZM
    F55 = np.diag([-betaf[0][0],-betaf[1][0],-betaf[2][0]]) #加速度計隨機模式矩陣
    F5 = np.hstack([F51, F52, F53, F54, F55])
    '''whole F'''
    F = np.vstack([F1, F2, F3, F4, F5])
    '''STM'''
    LO = np.eye(F.shape[0],F.shape[0]) + F*dt
    return LO

def KF_prediction (ERR_ST_old, LO, G, Q, P_old, dt):
    ERR_ST1 = LO[0:2,3:5].dot(ERR_ST_old[3:5])*(180/pi)
    ERR_ST2 = LO[2,5]*ERR_ST_old[5]
    ERR_ST3 = LO[3:6,6:9].dot(ERR_ST_old[6:9])*(pi/180) + LO[3:6,12:15].dot(ERR_ST_old[12:15])
    ERR_ST4 = LO[6:9,3:6].dot(ERR_ST_old[3:6])*(180/pi) + LO[6:9,9:12].dot(ERR_ST_old[9:12])
    ERR_ST5 = LO[9:12,9:12].dot(ERR_ST_old[9:12])
    ERR_ST6 = LO[12:15,12:15].dot(ERR_ST_old[12:15])
    ERR_ST_new = np.vstack([ERR_ST1, ERR_ST2, ERR_ST3, ERR_ST4, ERR_ST5, ERR_ST6])
    P_new = LO.dot(P_old).dot(LO.T) + (G*dt).dot(Q).dot((G*dt).T)
    return ERR_ST_new, P_new

def KF_correction (ERR_ST_new, P_new, Zk, Rk):
    Hk1 = np.hstack([np.eye(3,3), np.zeros([3,3]), np.zeros([3,9])])
    Hk2 = np.hstack([np.zeros([3,3]), np.eye(3,3), np.zeros([3,9])])
    Hk3 = np.hstack([np.zeros([2,3]),np.zeros([2,3]),np.eye(2,2),np.zeros([2,7])])
    Hk = np.vstack([Hk1,Hk2,Hk3])
    Kk = P_new.dot(Hk.T).dot(np.linalg.inv(Hk.dot(P_new).dot(Hk.T)+Rk))
    '''
    ERR_ST_new[0:2] = ERR_ST_new[0:2]*pi/180
    ERR_ST_new[6:9] = ERR_ST_new[6:9]*pi/180
    ERR_ST_new[9:12] = ERR_ST_new[9:12]*pi/180
    Zk[0:2] = Zk[0:2]*pi/180
    '''
    
    ERR_ST_new_c = ERR_ST_new + Kk.dot(Zk-Hk.dot(ERR_ST_new))
    P_new_c = (np.eye(Hk.shape[1],Hk.shape[1]) - Kk.dot(Hk)).dot(P_new) #tranditional form
    '''
    ERR_ST_new_c[0:2] = ERR_ST_new_c[0:2]*180/pi
    ERR_ST_new_c[6:9] = ERR_ST_new_c[6:9]*180/pi
    ERR_ST_new_c[9:12] = ERR_ST_new_c[9:12]*180/pi
    '''
    return ERR_ST_new_c, P_new_c