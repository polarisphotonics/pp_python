# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 22:17:51 2021

@author: User
"""

import numpy as np
from numpy import cos, sin, pi, sqrt, arctan, arcsin, arctan2
import random as rd
'''Direct cosine matrix'''
def DCM(p, r, y):
    '''
    Rz = np.array([
        [cos(y*pi/180), -sin(y*pi/180), 0],
        [sin(y*pi/180),  cos(y*pi/180), 0],
        [            0,              0, 1]])
    '''
    Rz = np.array([
        [cos(y*pi/180), sin(y*pi/180), 0],
        [-sin(y*pi/180),  cos(y*pi/180), 0],
        [            0,              0, 1]]) #y是AZ
    Rx = np.array([
        [1,             0,              0],
        [0, cos(p*pi/180), -sin(p*pi/180)],
        [0, sin(p*pi/180),  cos(p*pi/180)]])
    Ry = np.array([
        [ cos(r*pi/180), 0, sin(r*pi/180)],
        [             0, 1,             0],
        [-sin(r*pi/180), 0, cos(r*pi/180)]])
    R = Rz.dot(Rx).dot(Ry)
    RT = np.transpose(R)
    return R, RT

def ORI2Quaternion(p, r, y):
    [R,RT] = DCM(p,r,y)
    q4 = 0.5*sqrt(1+R[0][0]+R[1][1]+R[2][2])
    q1 = 0.25*(R[2][1]-R[1][2])/q4
    q2 = 0.25*(R[0][2]-R[2][0])/q4
    q3 = 0.25*(R[1][0]-R[0][1])/q4
    Q = np.array([[q1],[q2],[q3],[q4]])
    return Q, R

def Quaternion2DCM(Q):
    q1 = Q[0][0]
    q2 = Q[1][0]
    q3 = Q[2][0]
    q4 = Q[3][0]
    R = np.array([[q1**2-q2**2-q3**2+q4**2,  2*(q1*q2-q3*q4), 2*(q1*q3+q2*q4)],
                   [2*(q1*q2+q3*q4), -q1**2+q2**2-q3**2+q4**2, 2*(q2*q3-q1*q4)],
                   [2*(q1*q3-q2*q4), 2*(q2*q3+q1*q4), -q1**2-q2**2+q3**2+q4**2]])
    return R

def Quaternion_update(omg, Q_k, dt):
    theta_x = omg[0][0]*dt*pi/180
    theta_y = omg[1][0]*dt*pi/180
    theta_z = omg[2][0]*dt*pi/180
    theta = sqrt((theta_x)**2 + (theta_y)**2 + (theta_z)**2)\
        +0.000000000000000001
    S = (2/(theta))*sin(theta/2) ###########################改過theta單位喔
    C = 2*(cos(theta/2)-1)
    Mat = np.array([[C, S*theta_z, -S*theta_y, S*theta_x],
                    [-S*theta_z, C, S*theta_x, S*theta_y],
                    [S*theta_y, -S*theta_x, C, S*theta_z],
                    [-S*theta_x, -S*theta_y, -S*theta_z, C]])
    Q_k1 = Q_k + 0.5*Mat.dot(Q_k)
    R = Quaternion2DCM(Q_k1)
    ORI = DCM2Euler(R)
    return R, Q_k1, ORI

def DCM2Quaternion(R):
    q4 = 0.5*sqrt(1+R[0][0]+R[1][1]+R[2][2])
    q1 = 0.25*(R[2][1]-R[1][2])/q4
    q2 = 0.25*(R[0][2]-R[2][0])/q4
    q3 = 0.25*(R[1][0]-R[0][1])/q4
    Q = np.array([[q1],[q2],[q3],[q4]])
    return Q

def DCM2Euler(R):
    #Pitch = arctan(R[2][1]/sqrt((R[0][1])**2+(R[1][1])**2))*180/pi
    R = R + 0.0000000000000001*np.eye(3,3)*rd.random()
    Pitch = arcsin(R[2][1])*180/pi
    #Roll = -arctan(R[2][0]/R[2][2])*180/pi
    Roll = -arctan2(R[2][0],R[2][2])*180/pi
    Yaw = arctan2(R[0][1],R[1][1])*180/pi
    #Yaw = -arctan2(R[0][1],R[1][1])*180/pi
    '''
    if R[2][2] >0:
        Roll = Roll
    elif R[2][2] <0 and Roll <0:
        Roll = Roll #+ 180
    elif R[2][2] <0 and Roll >0:
        Roll = Roll #- 180

    if R[1][1] >0 and Yaw>0:
        Yaw = Yaw
    elif R[1][1] >0 and Yaw<0:
        Yaw = Yaw + 360
        #Yaw = Yaw + 180
    elif R[1][1] <0:
        Yaw = Yaw + 180
        #Yaw = Yaw + 360
    '''    
    ORI = np.array([[Pitch], [Roll], [Yaw]])
    return ORI