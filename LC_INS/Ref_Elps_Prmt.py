# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:48:43 2021

@author: User
"""
#import tool

from numpy import pi, sqrt

'''WGS84 parameter setting'''


#Major-axis 
Major_R = 6378137 #(unit: m)

#Flatness
Flat = 1/298.257223563 

#Eccentricity
Eccen = sqrt(Flat*(2-Flat))

#Semi-axis
Minor_R = Major_R*(1 - Flat) #(unit: m)

# rotation rate of the Earth
omg_iee = 7.292115*(10**(-5))*(180/pi) #(unit:deg/s)

#Parameters of Somigliana formula for determining ellisoidal gravity
SF1 = 9.7803267714
SF2 = 0.0052790414
SF3 = 0.0000232718
SF4 = -0.0000030876910891
SF5 = 0.0000000043977311
SF6 = 0.0000000000007211
