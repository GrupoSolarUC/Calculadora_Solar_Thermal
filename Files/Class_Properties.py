# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:08:53 2022

@author: Adrian Riebel
"""

import numpy as np
from numpy import polyfit, poly1d


class Properties:
    def __init__(self, glycol_percentage, P):
        self.glycol_percentage = glycol_percentage
        sat_data = np.load('Water_Glycol_Mixture_Saturation_Data.npz')
        sat_data = sat_data['G'+str(glycol_percentage)]
        for row in sat_data:
            if row[0] == P:
                self.T_sat = row[1]
                h_sat_liq = row[2]
                h_sat_vap = row[3]
                break
        h_vs_T_data = np.load('Water_Glycol_Mixture_Enthalpy_vs_Temp_Data.npz')
        h_vs_T_data = h_vs_T_data['G'+str(glycol_percentage)+'_P'+str(int(np.round(P*10,0)))]
        h_min = h_vs_T_data[0][1]
        h_to_add = 2*abs(h_min) + 50
        self.h_sat_liq = h_sat_liq + h_to_add
        self.h_sat_vap = h_sat_vap + h_to_add
        if glycol_percentage == 0:
            self.fluid_spec_heat = 4.184
            self.fluid_density = 1000
            self.fluid_thermal_conductivity = 2.263
            vapor_temps = []
            vapor_enthalpies = []
            for i in range(len(h_vs_T_data)):
                if h_vs_T_data[i][0] > self.T_sat:
                    vapor_temps.append(h_vs_T_data[i][0])
                    vapor_enthalpies.append(h_vs_T_data[i][1] + h_to_add)
            self.h_to_T_vap = poly1d(polyfit(vapor_enthalpies, vapor_temps,3))
            self.T_to_h_vap = poly1d(polyfit(vapor_temps, vapor_enthalpies,3))
        else:
            liquid_temps = []
            liquid_enthalpies = []
            for i in range(len(h_vs_T_data)):
                if h_vs_T_data[i][0] < self.T_sat:
                    liquid_temps.append(h_vs_T_data[i][0])
                    liquid_enthalpies.append(h_vs_T_data[i][1] + h_to_add)
            self.fluid_spec_heat = (liquid_enthalpies[len(liquid_enthalpies) - 1] - liquid_enthalpies[0])/(liquid_temps[len(liquid_temps) - 1] - liquid_temps[0])
            self.h_to_T_liq = poly1d(polyfit(liquid_enthalpies, liquid_temps,3))
            self.T_to_h_liq = poly1d(polyfit(liquid_temps, liquid_enthalpies,3))
    def h_to_T(self, h):
        if self.glycol_percentage == 0:
            if h > self.h_sat_vap:
                return self.h_to_T_vap(h)
            if h < self.h_sat_liq:
                return self.T_sat - (self.h_sat_liq - h)/self.fluid_spec_heat
            return self.T_sat
        else:
            if h >= self.h_sat_liq:
                return self.T_sat
            else:
                return self.h_to_T_liq(h)
    def T_to_h(self, T):
        if self.glycol_percentage == 0:
            if T > self.T_sat:
                return self.T_to_h_vap(T)
            elif T == self.T_sat:
                raise ValueError('Specific enthalpy cannot be computed at saturation temperature: '
                                 +str(self.T_sat)+' °C')
            else:
                return self.h_sat_liq - self.fluid_spec_heat*(self.T_sat - T)
        else:
            if T > self.T_sat:
                raise ValueError('Temperature surpasses the saturation temperature ( = '+str(self.T_sat)+' ). This calculation can only be done with pure water')
            elif T == self.T_sat:
                raise ValueError('Specific enthalpy cannot be computed at saturation temperature: '
                                 +str(self.T_sat)+' °C')
            else:
                return self.T_to_h_liq(T)

    def quality(self, h):
        if h > self.h_sat_vap:
            return 1
        if h < self.h_sat_liq:
            return 0
        return (h - self.h_sat_liq)/(self.h_sat_vap - self.h_sat_liq)

