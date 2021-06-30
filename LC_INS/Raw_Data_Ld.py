# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 15:11:44 2021

@author: User
"""


#import tool
import numpy as np
'''
Illustraction for data format:
--------------------------------------------------------------------------------------------------------------------------------------------------
#2021-04-22 12:13:46 (start time)															
#dt	 SRS200_wz	 PP_wz	 Nano33_wx	 Nano33_wy	 Nano33_wz	 Adxl355_ax  Adxl355_ay	 Adxl355_az	 speed	 VBOX	 data_T	lat	lon	sat_N	GPS_h Velo_u
#s	       DPS	   DPS	       DPS	       DPS	       DPS         	  g           g	          g	   m/s	  m/s	   code	deg	deg	  num	    m   m/s
.            .       .           .           .           .            .           .           .     .      .          .  .   .      .       .    .
.            .       .           .           .           .            .           .           .     .      .          .  .   .      .       .    .
.            .       .           .           .           .            .           .           .     .      .          .  .   .      .       .    .
#2021-04-22 12:31:05(end time)
--------------------------------------------------------------------------------------------------------------------------------------------------
dt: 相對時間起點的時間間隔，後續需化做UTC+8的時間，並以秒為單位(e.g. dt+12*3600+13*60+46)
SRS200_wz: 為戰術應用中階FOG所獲Z軸角速度，單位degree per sec (DPS)
PP_wz: 為PolarisPhotonics(PP)開發之FOG所獲Z軸角速度，單位degree per sec (DPS)
Nano33_wx,Nano33_wy,Nano33_wz: 為MEMS Nano-33型號IMU所獲X,Y,Z軸角速度，單位degree per sec (DPS)
Adxl355_ax, Adxl355_ay, Adxl355_az: 為MEMS Adxl355型號accelerometer所獲X,Y,Z軸加速度，單位g(後續還須根據所在地球位置換算成m/s/s)
speed: unknown
VBOX: VBOX輸出之GNSS velocity
data_T: unknown
lat: VBOX輸出之GNSS緯度
lon: VBOX輸出之GNSS經度
sat_N: VBOX輸出之衛星顆數(GPS+Glonass)
GPS_h: VBOX輸出之GNSS高度(橢球高)
Velo_u: VBOX輸出之高程速度
'''
NAME = 'RT4_OBS'
Raw_data = np.loadtxt(NAME+'.txt', comments='#', delimiter=',')
num_pen = Raw_data.shape[0]
#num_pen = 20000