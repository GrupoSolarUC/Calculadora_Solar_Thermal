# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:08:53 2022

@author: Adrian
"""

from math import pi, exp


class Storage_Tank:
    '''
    One-dimensional model for a vertical, cylindrical, stratified tank for thermal energy storage.
    This model divides the tank into nodes which are vertically arranged. Within each node, the temperature is assumed to be uniform.
    The model assumes that the fluid stored has properties (density, specific heat capacity, thermal conductivity) which are undependent from its temperature.
    The model allows two inlet-outlet pairs. The flow leaving the tank through an outlet is automatically equal to the flow entering the tank through its associated inlet, so that the mass and volume in the tank are maintained.
    
    Initialization parameters:
        - Volume: Total volume of the tank (m^3).
        - AR: Aspect ratio ( i.e. height/diameter ) (no units)
        - top_loss_coeff: thermal loss coefficient in the top surface of the tank ( kJ / ( hr * m^2 *K ) ).
        - edge_loss_coeff: thermal loss coefficient in the lateral surface of the tank ( kJ / ( hr * m^2 *K ) ).
        - bottom_loss_coeff: thermal loss coefficient in the bottom surface of the tank ( kJ / ( hr * m^2 *K ) ).
        - N_nodes: number of nodes into which the tank is divided (no units).
        - inlet_1_node: node where the first inlet to the tank is located (no units).
        - outlet_1_node: node where the first outlet to the tank (associated to inlet 1) is located (no units).
        - inlet_2_node: node where the second inlet to the tank is located (no units).
        - outlet_2_node: node where the second outlet to the tank (associated to inlet 2) is located (no units).
        - fluid_props: object of the class "Properties_water" to define the properties of the fluid.
        - initial_temperatures: dictionary with the initial temperatures of each node. It has the format: {1: T1, 2: T2, 3: T3, ..... , N: TN } where N is the number of nodes (node 1 is the top node) and T1, T2, etc are the temperatures of the nodes in °C.
        - tolerance: maximum discrepancy between successive iterations of the tank which is allowed to yield the tank outputs (no units).
        - max_iterations: maximum number of iterations permitted until the outputs of the tank are yielded despite no convergence has been reached (no units).
        
    Attributes (See Initialization parameters for the ones that have no description):
        - self.cross_Area: cross-sectional area of the tank
        - self.N_nodes
        - self.node_temperatures: dictionary contanining the current temperatures. It is initialized as "initial_temperatures" (See Initialization parameters)
        - self.top_loss_coeff
        - self.edge_loss_coeff
        - self.bottom_loss_coeff
        - self.node_height: height of each node
        - self.nodal_edge_area: lateral area of each node
        - self.node_volume: volume of each node
        - self.inlets_and_outlets: dictionary to store the nodes where the inlets and outlets are located. It has the format: {'inlet_1': inlet_1_node, 'outlet_1': outlet_1_node, 'inlet_2': inlet_2_node, 'outlet_2': outlet_2_node} (See initialization parameters)
        - self.max_iterations
        - self.fluid_props
        - self.tolerance
        - self.new_temperatures: dictionary used to store the last computed outputs, in order to update the temperatures of the tank once convergence is achieved in the rest of the system. This attribute is not defined until the self.compute_outputs method is used for the first time (See that method for more detail)
    
    Methods:
        - self.convergence: evaluates whether convergence has been achieved in order to yield the solutions.
        - self.mix_nodes: mixes unstable nodes at the end of the time step
        - self.iteration: executes an iteration to find the solution of the current time step
        - self.compute_outputs: executes several iterations until convergence is achieved. Then it returns the resulting specific enthalpies at the outlets.
        - self.update_temperature: updates the temperatures of the tank to the ones determined in the last output-computation
    '''
    def __init__(self, volume, AR, top_loss_coeff, edge_loss_coeff, bottom_loss_coeff, N_nodes,
                 inlet_1_node, outlet_1_node, inlet_2_node, outlet_2_node, fluid_props,
                 initial_temperatures, tolerance = 1e-8, max_iterations = 100):
        diameter = (4*volume/(pi*AR))**(1/3)
        height = diameter*AR
        radius = diameter/2
        self.cross_Area = pi*radius**2
        self.N_nodes = N_nodes ## number of nodes
        self.node_temperatures = initial_temperatures
        self.top_loss_coeff = top_loss_coeff ## loss coefficient at the top
        self.edge_loss_coeff = edge_loss_coeff ## loss coefficient at the edge
        self.bottom_loss_coeff = bottom_loss_coeff ## loss coefficient at the bottom
        self.node_height = height/N_nodes ## height of each node
        self.nodal_edge_area = 2*pi*radius*self.node_height ## lateral area of each node
        self.node_volume = volume/N_nodes ## volume of each node
        ## self.inlets_and_outlets: dictionary to register the nodes where the inlets and outlets are
        self.inlets_and_outlets = {'inlet_1': inlet_1_node, 'outlet_1': outlet_1_node,
                                   'inlet_2': inlet_2_node, 'outlet_2': outlet_2_node}
        self.max_iterations = max_iterations
        self.fluid_props = fluid_props
        self.tolerance = tolerance
    ## convergence: Function that tests the convergence considering successive
    ## iterations of the tank.
    ## Parameters:
    ##      - dict1: dictionary with temperatures of the previous iteration
    ##      - dict2: dictionary with temperatures of the current iteration
    ##      - tolerance: number to define how small the difference between the iterations can be
    def convergence(self, dict1, dict2):
        '''
        Function that evaluates the results of two successive iterations and returns a boolean value depending on whether the iterations have converged.
        
        Parameters:
            - dict1 : dictionary that contains the results of the previous iteration (mean and final temperatures in °C).
            - dict2 : dictionary that contains the results of the current iteration (mean and final temperatures in °C).
        
        Returns:
            - Boolean value which is True if convergence has been achieved and False otherwise.
        '''
        assert len(dict1) == len(dict2)
        N_nodes = len(dict1)
        List1 = [ abs((dict1[i]['T_ave'] -
                       dict2[i]['T_ave'])/(dict1[i]['T_ave']+273.15)) < self.tolerance
                 for i in range(1, N_nodes+1) ]
        List2 = [ abs((dict1[i]['T_final'] -
                       dict2[i]['T_final'])/(dict1[i]['T_final']+273.15)) < self.tolerance
                 for i in range(1, N_nodes+1) ]
        return all(List1) and all(List2)
    ## self.mix_nodes: function that mixes the nodes that are unstable. I.e. when some node is hotter than the node above it,
    ## this function mixes those nodes. The process is repeated until there is no instabilities.
    def mix_nodes(self, dict1):
        '''
        Function that mixes the nodes at the end of the time step if there are any unstabilities in the tank (i.e. a node is hotter than the node above it)
        
        Parameters:
            A dictionary returned by the "iteration" method, with the format: {1: {'T_ave': T_ave_1, 'T_final': T_final_1}, 2: {'T_ave': T_ave_2, 'T_final': T_final_2}, ... , N: {'T_ave': T_ave_N, 'T_final': T_final_N} } (See the "return" description of the "iteration" method)
        
        Returns:
            A dictionary with the same format as the input, but each value of 'T_final' has been updated as described. The values of 'T_ave' are not modified.
        '''
        mixed_nodes = [ {'original_node_numbers': [i], 'T_final': dict1[i]['T_final'],
                         'mass':1} for i in dict1 ]
        while True:
            List1 = [ mixed_nodes[i]['T_final'] >= mixed_nodes[i+1]['T_final'] for i in range(len(mixed_nodes) - 1) ]
            if all(List1):
                break
            for i in range(len(mixed_nodes) - 1):
                if mixed_nodes[i]['T_final'] < mixed_nodes[i+1]['T_final']:
                    total_mass = mixed_nodes[i]['mass'] + mixed_nodes[i+1]['mass']
                    T_final = ((mixed_nodes[i]['mass']/total_mass)*mixed_nodes[i]['T_final'] +
                               (mixed_nodes[i+1]['mass']/total_mass)*mixed_nodes[i+1]['T_final'])
                    nodes_list = mixed_nodes[i]['original_node_numbers'] + mixed_nodes[i+1]['original_node_numbers']
                    mixed_nodes[i] = {'original_node_numbers': nodes_list, 'T_final': T_final,
                                      'mass': total_mass}
                    mixed_nodes.pop(i+1)
                    break
        R1 = {}
        for i in range(len(mixed_nodes)):
            for node_number in mixed_nodes[i]['original_node_numbers']:
                R1[node_number] = {'T_ave': dict1[node_number]['T_ave'], 'T_final': mixed_nodes[i]['T_final']}
        return R1
    ## self.iteration: it executes a single iteration to determine the temperatures of the nodes at the end of
    ## the time step and the average temperatures during the time step
    ## It usually must be called more than one time to converge to definitve outputs given a certain input.
    def iteration(self, inlet_1_flow_rate, inlet_1_T, inlet_2_flow_rate, inlet_2_T,
                  T_amb, previous_result, fluid_spec_heat, fluid_thermal_conductivity,
                  node_thermal_capacity, time_step):
        '''
        Function that executes one iteration to determine the mean temperatures of the nodes during the time step and the temperatures of the nodes at the end of the time step.

        Parameters:
            - inlet_1_flow_rate: mass flow through inlet 1 (kg/hr)
            - inlet_1_T: temperature of the flow entering through inlet 1 (°C)
            - inlet_2_flow_rate: mass flow through inlet 2 (kg/hr)
            - inlet_2_T: temperature of the flow entering through inlet 2 (°C)
            - T_amb: ambient temperature outside the tank (°C)
            - previous_result: temperatures determined in the previous iteration (°C). This dictionary has the same format of the dictionary that this function returns.
            - fluid_spec_heat: specific heat capacity of the fluid ( kJ / (kg * K) )
            - fluid_thermal_conductivity: thermal conductivity of the fluid ( kJ / (hr * m * K) )
            - node_thermal_capacity: heat capacity of the node ( kJ / K )
            - time_step: length of the time step that is being used (hr)

        Returns
            A dictionary with the following format: {1: {'T_ave': T_ave_1, 'T_final': T_final_1}, 2: {'T_ave': T_ave_2, 'T_final': T_final_2}, ... , N: {'T_ave': T_ave_N, 'T_final': T_final_N} }
            where N is the number of nodes, T_ave_i is the average temperature of the i-th node during the time step, and T_final_i is the temperature of that node at the end of the time step 
        '''
        Sol = {}
        for i in range(1, self.N_nodes + 1):
            ## values of a and b are computed, according to the mathematical reference of type 158 of TRNSYS
            a = -self.nodal_edge_area*self.edge_loss_coeff
            b = self.nodal_edge_area*self.edge_loss_coeff*T_amb
            ## if the node has inlets or outlets, flow_in and flow_out are updated
            ## this also modifies the values of a and b
            flow_in  = 0
            if self.inlets_and_outlets['inlet_1'] == i:
                flow_in = flow_in + inlet_1_flow_rate
                a = a - inlet_1_flow_rate*fluid_spec_heat
                b = b + inlet_1_flow_rate*fluid_spec_heat*inlet_1_T
            if self.inlets_and_outlets['inlet_2'] == i:
                flow_in = flow_in + inlet_2_flow_rate
                a = a - inlet_2_flow_rate*fluid_spec_heat
                b = b + inlet_2_flow_rate*fluid_spec_heat*inlet_2_T
            flow_out = 0 ## flow going out of the node through outlets of the tank
            if self.inlets_and_outlets['outlet_1'] == i:
                flow_out = flow_out + inlet_1_flow_rate
            if self.inlets_and_outlets['outlet_2'] == i:
                flow_out = flow_out + inlet_2_flow_rate
            ## Case 1: Top node of the tank
            if i == 1:
                #flow_below = flow coming from the node below (if positive) or going to the node below (if negative)
                flow_below = flow_out - flow_in
                a = a - self.cross_Area*self.top_loss_coeff
                a = a - fluid_thermal_conductivity*self.cross_Area/self.node_height
                b = b + self.cross_Area*self.top_loss_coeff*T_amb
                b = b + fluid_thermal_conductivity*self.cross_Area*previous_result[i+1]['T_ave']/self.node_height
                if flow_below > 0:
                    a = a - flow_below*fluid_spec_heat
                    b = b + flow_below*fluid_spec_heat*previous_result[i+1]['T_ave']
                ## flow_above = -flow_below for the next node
                flow_above = -flow_below
            ## Case 2: Bottom node of the tank
            elif i == self.N_nodes:
                a = a - self.cross_Area*self.bottom_loss_coeff
                a = a - fluid_thermal_conductivity*self.cross_Area/self.node_height
                b = b + self.cross_Area*self.bottom_loss_coeff*T_amb
                b = b + fluid_thermal_conductivity*self.cross_Area*previous_result[i-1]['T_ave']/self.node_height
                ## since this is the bottom node, a warning is printed if 
                ## the flow coming from or going to a node below is larger than zero
                if abs(flow_above + flow_in - flow_out) > 0.0001:
                    print('Warning: A problem has occurred in the internal flow balance of the tank')
                if flow_above > 0:
                    a = a - flow_above*fluid_spec_heat
                    b = b + flow_above*fluid_spec_heat*previous_result[i-1]['T_ave']
            ## Case 3: Intermediate nodes of the tank
            else:
                ## flow_below corresponds to the flow coming from the node below (if positive)
                ## or going to the node below (if negative)
                ## flow_above has been defined when computing the previous node (node above)
                flow_below = flow_out - flow_in - flow_above
                a = a - 2*fluid_thermal_conductivity*self.cross_Area/self.node_height
                b = b + fluid_thermal_conductivity*self.cross_Area*previous_result[i+1]['T_ave']/self.node_height
                b = b + fluid_thermal_conductivity*self.cross_Area*previous_result[i-1]['T_ave']/self.node_height
                if flow_above > 0:
                    a = a - flow_above*fluid_spec_heat
                    b = b + flow_above*fluid_spec_heat*previous_result[i-1]['T_ave']
                if flow_below > 0:
                    a = a - flow_below*fluid_spec_heat
                    b = b + flow_below*fluid_spec_heat*previous_result[i+1]['T_ave']
                ## flow_above = -flow_below for the next node
                flow_above = -flow_below
            a = a/node_thermal_capacity
            b = b/node_thermal_capacity
            ## T_final: temperature of the node at the end of the time_step
            ## T_ave: average temperature during the time_step.
            T_final = (self.node_temperatures[i]+b/a)*exp(a*time_step)-b/a
            T_ave = (1/(a*time_step))*(self.node_temperatures[i]+b/a)*(exp(a*time_step)-1)-b/a
            Sol[i] = {'T_ave': T_ave, 'T_final': T_final}
        return Sol
    ## self.compute_outputs: given the inputs of the tank, this function calls the "iteration" function
    ## until convergence is achieved.
    def compute_outputs(self, inlet_1_flow_rate, inlet_1_h, inlet_2_flow_rate, inlet_2_h,
                        T_amb, time_step):
        '''
        Function that calls the "iteration" function several times until the temperatures of the nodes reach convergence.
        
        This function computes the outputs of the tank but it does not update the tank temperatures, because the rest of the system must achieve convergence as well before the tank's temperatures are updated.

        Parameters:
            - inlet_1_flow_rate: mass flow through inlet 1 (kg/hr)
            - inlet_1_h: specific enthalpy of the flow entering through inlet 1 (°C)
            - inlet_2_flow_rate: mass flow through inlet 2 (kg/hr)
            - inlet_2_h: specific enthalpy of the flow entering through inlet 2 (°C)
            - T_amb: ambient temperature outside the tank (°C)
            - time_step: length of the time step that is being used (hr)

        Returns:
            a dictionary  with the following keys:
                - 'outlet_1_h': specific enthalpy of the flow leaving the tank through outlet 1
                - 'outlet_2_h': specific enthalpy of the flow leaving the tank through outlet 2
            
            This function also defines a new attribute: self.new_temperatures which is a dictionary to store the temperatures computed, in order to update the self.temperatures attribute of the model once convergence is achieved in the rest of the system. It has the format: {1: T_final_1, 2: T_final_2, .... , N: T_final_N} where T_final_i is the temperature of the i-th node at the end of the time step.
        '''
        ## the tempratures to initialize the iterations are the temperatures of the nodes at the beginning of the time step
        previous_result = {i: {'T_ave':self.node_temperatures[i], 'T_final': self.node_temperatures[i]}
                           for i in range(1, self.N_nodes + 1)}
        fluid_spec_heat = self.fluid_props.fluid_spec_heat
        fluid_density = self.fluid_props.fluid_density
        fluid_thermal_conductivity = self.fluid_props.fluid_thermal_conductivity
        node_thermal_capacity = self.node_volume*fluid_density*fluid_spec_heat ## thermal capacity of each node
        inlet_1_T = self.fluid_props.h_to_T(inlet_1_h)
        inlet_2_T = self.fluid_props.h_to_T(inlet_2_h)
        it = 0
        ## loop to call the iteration function
        while True:
            new_result = self.iteration(inlet_1_flow_rate, inlet_1_T, inlet_2_flow_rate, inlet_2_T, T_amb,
                                        previous_result, fluid_spec_heat, fluid_thermal_conductivity,
                                        node_thermal_capacity, time_step)
            it = it + 1
            ## when convergence is achieved, the loop is ended
            if self.convergence(previous_result, new_result):
                maxIt = False
                break
            if it == self.max_iterations:
                maxIt = True
                break
            previous_result = new_result
        ## If some node is hotter than the node above it, the nodes are mixed until this problem is solved.
        new_result = self.mix_nodes(new_result)
        ## self.new_temperatures is a provisional dictionary to save the temperatures that the algorithm determined.
        ## This temperatures are set as definitive temperatures of the tank only if the rest of the system has converged as well.
        ## If the inputs to the tank change during the current time step because no convergence has been achieved,
        ## the outputs must be computed again beginning with the temperatures stored in self.node_temperatures
        self.new_temperatures = {i: new_result[i]['T_final'] for i in new_result}
        ## extract outputs from the resulting temperatures
        outlet_1_T = new_result[self.inlets_and_outlets['outlet_1']]['T_ave']
        outlet_2_T = new_result[self.inlets_and_outlets['outlet_2']]['T_ave']
        heat_rate_to_demand = inlet_2_flow_rate*fluid_spec_heat*(outlet_2_T - inlet_2_T)
        Heat_Loss_Rate = sum([ self.nodal_edge_area*self.edge_loss_coeff*(new_result[i]['T_ave'] - T_amb)
                              for i in new_result])
        Heat_Loss_Rate = Heat_Loss_Rate + self.cross_Area*self.top_loss_coeff*(new_result[1]['T_ave'] - T_amb)
        Heat_Loss_Rate = Heat_Loss_Rate + self.cross_Area*self.bottom_loss_coeff*(new_result[self.N_nodes]['T_ave'] - T_amb)
        return {'outlet_1_h': self.fluid_props.T_to_h(outlet_1_T),
                'outlet_2_h': self.fluid_props.T_to_h(outlet_2_T),
                'Q_demand': heat_rate_to_demand, 'Q_loss': Heat_Loss_Rate, 'maxIt': maxIt }
    ## self.update_temperature: only when confirming that the inputs and ouputs of the whole system have converged,
    ## the temperature of the tank is updated to the value that it previously determined
    def update_temperature(self):
        '''
        Function to update the temperature of the tank according to the last outputs computed, when convergence has been achieved in the whole system.
        
        No parameters, no returns.
        '''
        self.node_temperatures = self.new_temperatures
