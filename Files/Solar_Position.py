# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 23:31:26 2022

@author: adria
"""

from scipy.interpolate import interp1d

def zenith_function(zenith_list):
    assert len(zenith_list)%8760 == 0
    factor = int(len(zenith_list)/8760)
    time = [ t/factor for t in range(len(zenith_list)) ]
    time.append(8760)
    zenith_list.append(zenith_list[0])
    return interp1d(time, zenith_list)

def azimuth_function(azimuth_list):
    assert len(azimuth_list)%8760 == 0
    factor = int(len(azimuth_list)/8760)
    time = [ t/factor for t in range(len(azimuth_list)) ]
    i = 0
    while i <  len(azimuth_list) - 1:
        if azimuth_list[i + 1] - azimuth_list[i] > 200:
            if azimuth_list[i] == 0 and azimuth_list[i + 1] == 360:
                i = i + 1
            else:
                delta_t = time[i + 1] - time[i]
                delta_t_1 = azimuth_list[i]*delta_t/(azimuth_list[i] + (360 - azimuth_list[i + 1]))
                time.insert(i + 1, time[i] + delta_t_1)
                azimuth_list.insert(i + 1, 0)
                time.insert(i + 2, time[i] + delta_t_1)
                azimuth_list.insert(i + 2, 360)
                i = i + 3
        elif azimuth_list[i] - azimuth_list[i + 1] > 200:
            if azimuth_list[i] == 360 and azimuth_list[i + 1] == 0:
                i = i + 1
            else:
                delta_t = time[i + 1] - time[i]
                delta_t_1 = (360 - azimuth_list[i])*delta_t/((360 - azimuth_list[i]) + azimuth_list[i + 1])
                time.insert(i + 1, time[i] + delta_t_1)
                azimuth_list.insert(i + 1, 360)
                time.insert(i + 2, time[i] + delta_t_1)
                azimuth_list.insert(i + 2, 0)
                i = i + 3
        else:
            i = i + 1
    time.append(8760)
    azimuth_list.append(azimuth_list[0])
    return interp1d(time, azimuth_list)