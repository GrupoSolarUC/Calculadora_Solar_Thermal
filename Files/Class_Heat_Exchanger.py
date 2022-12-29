# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:08:53 2022

@author: Adrian
"""

class Heat_Exchanger:
    '''
    Model for a constant effectiveness heat exchanger. The model considers energy losses due to vapor generation only at the load-side of the heat exchanger.
    When vapor is produced due to boiling in the load-side of the heat exchanger, a flow of mains water is introduced in the inlet of the load-side to compensate the vapor generated.
    
    Initialization parameters/Attributes:
        - eff: effectiveness of the heat exchanger
        - fluid_props_source: A "Properties" class object to define the thermodynamic properties of the energy source fluid.
        - fluid_props_load: A "Properties" class object to define the thermodynamic properties of the load fluid.
        - tolerance: The maximum discrepancy between successive iterations that is allowed in order to yield the outputs of the model.
        - maxIt: Maximum numer of iterations allowed until the outputs of the heat exchanger are yielded despite no convergence has been achieved.
    
    Methods:
        - self.convergence: evaluates whether convergence has been achieved.
        - self.compute_outputs: computes the outlet enthalpies of both flows and other relevant results.
    
    Exceptions:
        - The model yields a ValueError when the temperature of the fluid entering the load-side of the collector is higher than the saturation temperature of the source fluid (this error does not happen viseversa, in which case the vapor generation may need to be computed).
    '''
    def __init__(self, eff, fluid_props_source, fluid_props_load, tolerance = 1e-6, maxIt = 50):
        self.eff = eff
        self.fluid_props_source = fluid_props_source
        self.fluid_props_load = fluid_props_load
        self.tolerance = tolerance
        self.maxIt = maxIt
    def convergence(self, m_mains, m_vap):
        '''
        Function that evaluates the convergence of the iterations (iterations are only necessary when vapor is being generated in the load side of the heat exchanger).
        In each iteration, a flow of mains water is introduced into the load-side inlet to compensate the mass loss due to vapor expulsion.
        In order to achieve convergence, the mass flow of water being introduced and the mass flow of vapor being produced must have a discrepancy smaller than the "tolerance" attribute.

        Parameters:
            - m_mains: mass flow of mains water that was introduced into the collector at the beginning of the iteration (kg/hr).
            - m_vap: mass flow of vapor produced due to the heat of the collectors (kg/hr).

        Returns:
            A boolean value which is True if convergence has been achieved, and False otherwise.
        '''
        if m_mains == m_vap:
            return True
        m_max = max([ m_mains, m_vap ])
        return abs(m_mains - m_vap)/m_max < self.tolerance
    def compute_outputs(self, m_in_source, h_in_source, m_in_load, h_in_load, h_mains):
        '''
        Functions that computes and yields the outputs of the heat exchanger, given the input flows.

        Parameters:
            - m_in_source: mass flow entering the "source" side of the heat exchanger (kg/hr)
            - h_in_source: specific enthalpy of the "source" fluid at the inlet (kJ/kg)
            - m_in_load: mass flow entering the "load" side of the heat exchanger (kg/hr)
            - h_in_load: specific enthalpy of the fluid entering the load-side (kJ/kg)
            - h_mains: enthalpy of the mains water flow that would be added to the load side if there is vapor production. (kJ/kg)

        Raises:
            - ValueError: It happens when the temperature of the fluid entering the load-side of the collector is higher than the saturation temperature of the source fluid (this error does not happen viseversa, in which case the vapor generation may need to be computed).

        Returns
            A dictionary with the keys:
                - 'h_out_load': specific enthalpy of the flow leaving the source side of the heat exchanger. (kJ/kg)
                - 'h_out_source': specific enthalpy of the flow leaving the load side of the heat exchanger. (kJ/kg)
                - 'm_vap': flow mass of vapor being released. (kg/hr)
                - 'Q_useful': all heat that is not used to produce vapor, and is thus useful. (kJ/hr)
                - 'Q_waste': heat that is used to produce vapor. (kJ/hr)
                - 'iterations': number of iterations that were necessary in order to achieve convergence. (no units)
                - 'maxIt': bolean value which is True if the maximum number of allowed iterations was reached without achieving convergence, and False otherwise.
        '''
        if m_in_source == 0 or m_in_load == 0:
            return {'h_out_load': h_in_load,
                    'h_out_source': h_in_source,
                    'm_vap': 0,
                    'Q_useful': 0,
                    'Q_waste': 0,
                    'iterations': 0,
                    'maxIt': False}
        T_in_source = self.fluid_props_source.h_to_T(h_in_source)
        m_mains = 0
        it = 0
        while True:
            new_m_in_load = m_in_load + m_mains
            new_h_in_load = (m_in_load*h_in_load + m_mains*h_mains)/new_m_in_load
            new_T_in_load = self.fluid_props_load.h_to_T(new_h_in_load)
            if new_T_in_load >= self.fluid_props_source.T_sat:
                raise ValueError('Load temperature entering heat exchanger is higher than the saturation temperature of the source fluid')
            h_out_source_hyp = self.fluid_props_source.T_to_h(new_T_in_load)
            if T_in_source >= self.fluid_props_load.T_sat:
                h_out_load_hyp = self.fluid_props_load.h_sat_liq + self.fluid_props_load.fluid_spec_heat*(T_in_source - self.fluid_props_load.T_sat)
            else:
                h_out_load_hyp = self.fluid_props_load.T_to_h(T_in_source)
            if T_in_source > new_T_in_load:
                Q_max = min([ m_in_source*(h_in_source - h_out_source_hyp),
                             new_m_in_load*(h_out_load_hyp - new_h_in_load) ])
            else:
                Q_max = -min([ m_in_source*(h_out_source_hyp - h_in_source),
                              new_m_in_load*(new_h_in_load - h_out_load_hyp) ])
            h_out_source = h_in_source - self.eff*Q_max/m_in_source
            h_out_load = new_h_in_load + self.eff*Q_max/new_m_in_load
            m_vap = new_m_in_load*self.fluid_props_load.quality(h_out_load)
            it = it + 1
            if self.convergence(m_mains, m_vap):
                maxIt = False
                break
            if it == self.maxIt:
                maxIt = True
                break
            Q_useful = m_in_load*(self.fluid_props_load.h_sat_liq - h_in_load)
            Q_actual = new_m_in_load*(h_out_load - new_h_in_load)
            Q_waste = Q_actual - Q_useful
            m_mains = Q_waste/(self.fluid_props_load.h_sat_vap - h_mains)
        if m_vap > 0:
            return {'h_out_load': self.fluid_props_load.h_sat_liq,
                    'h_out_source': h_out_source,
                    'm_vap': m_vap,
                    'Q_useful': Q_useful,
                    'Q_waste': Q_waste,
                    'iterations': it,
                    'maxIt': maxIt}
        else:
            return {'h_out_load': h_out_load,
                    'h_out_source': h_out_source,
                    'm_vap': 0,
                    'Q_useful': new_m_in_load*(h_out_load - new_h_in_load),
                    'Q_waste': 0,
                    'iterations': it,
                    'maxIt': maxIt}