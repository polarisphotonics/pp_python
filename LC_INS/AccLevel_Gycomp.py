# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 17:14:06 2021

@author: User
"""
#import tool

from numpy import sin, cos, arctan, pi, sqrt

def Acc_Leveling(AVG_fx, AVG_fy, AVG_fz, g_start):
    Pitch0 = arctan(AVG_fy/sqrt(AVG_fx**2 + (AVG_fz + g_start)**2))*(180/pi)
    Roll0 = arctan(-AVG_fx/(AVG_fz + g_start))*(180/pi)
    return Pitch0, Roll0
    

def Gyrocompassing(AVG_wx, AVG_wy, AVG_wz, Pitch0, Roll0, Lat0):
    omg_iee = 7.292115*(10**(-5))*(180/pi) #(unit:deg/s)
    AVG_wy = AVG_wy+omg_iee*cos(Lat0*pi/180)
    AVG_wz = AVG_wz+omg_iee*sin(Lat0*pi/180)
    Yaw0 = arctan((AVG_wx*cos(Roll0*pi/180)+AVG_wz*sin(Roll0*pi/180))\
                 /(AVG_wx*sin(Pitch0*pi/180)*cos(Roll0*pi/180)+\
                   AVG_wy*cos(Pitch0*pi/180)-AVG_wz*cos(Roll0*pi/180)\
                       *sin(Pitch0*pi/180)))*(180/pi)
    return Yaw0