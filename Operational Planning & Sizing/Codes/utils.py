# This file contains utility functions. You mostly do not want to touch this file.

from datetime import timedelta, datetime
from pyomo.opt import SolverFactory
from pyomo.opt import ProblemFormat
import pandas as pd
from parameters import RESOLUTION_IN_MINUTES, N_PERIODS,FUEL_CO2,GRID_CO2,PV_CO2,STORAGE_CO2,GENSET_CO2,INVERTER_CO2, INVESTMENT_HORIZON 
from parameters import EV_CAPACITY, EV_TARGET_SOC, BONUS, GRID_CONNECTED#, SIZING
import matplotlib.pyplot as plt
import numpy as np


from pyomo.core.kernel import value

class SizingData():
    def __init__(self, start_date=datetime(2021, 1, 1, 0, 0, 0), scen=0, yearly_cons=0):
        """
        :param start_date: The start date you want to consider for sizing.
        """
        # Parsing function
        def date_parse(x): 
            return datetime.strptime(x, '%Y-%m-%d %H:%M:%S')  
        
        # Columns are named "C" (consumption), "P" (production) and "EV" (electric vehicle)
        data_df = pd.read_csv("profile.csv", index_col="DateTime",
                              parse_dates=True, date_format='%Y-%m-%d %H:%M:%S')
        self.start_date = start_date
        self.data_df = data_df
        self.yearly_cons = yearly_cons
        
    def consumption(self, p):
        """
        :param p: period index, starts at 1
        :return: instantaneous power consumption [W]
        """
        datetime = self.start_date + \
            timedelta(hours=(p-1) * RESOLUTION_IN_MINUTES / 60)
        # DEBUG print(datetime.strftime('%Y-%m-%d %H:%M:%S'))
        return self.data_df.loc[datetime]["load"] * self.yearly_cons * 60 / RESOLUTION_IN_MINUTES

    def PV_generation(self, p):
        """
        :param p: period index, starts at 1
        :return: PV generation per unit of capacity [W]
        """
        datetime = self.start_date + \
            timedelta(hours=(p-1) * RESOLUTION_IN_MINUTES / 60)
        return self.data_df.loc[datetime]["PV"].clip(min=0)
    
    def EV_connected(self,p):
        #datetime = self.data_df.index[0] + \
        datetime = self.start_date + \
            timedelta(hours=(p-1) * RESOLUTION_IN_MINUTES / 60)
        return 1 if self.data_df.loc[datetime]["EV"].clip(min=0) > 0 else 0
    
    def EV_profile(self,p):
        datetime = self.data_df.index[0] + \
            timedelta(hours=(p-1) * RESOLUTION_IN_MINUTES / 60)
        return self.data_df.loc[datetime]["EV"].clip(min=0)


def extract_EV_data(data, yearly_km=0):
    n_connection = 0 
    previous = False
    #EV_connected = np.zeros(N_PERIODS)
    EV_initial_SOC = []
    t_arr = []
    t_dep = []
    for t in np.arange(1,N_PERIODS+1):
        if previous == False and data.EV_profile(t) > 0:
            n_connection += 1                                                   #Increase the number of connection
            k = yearly_km/5e3                                                   # Default profile is for 5k km 
            EV_initial_SOC.append(EV_TARGET_SOC-data.EV_profile(t)*k)           # EV_initial_SOC = target-energy_needed
                
            t_arr.append(t)
            previous = True
            #EV_connected[t-1] = 1
        elif data.EV_profile(t) > 0 and previous == True:
            #EV_connected[t-1] = 1
            if t == N_PERIODS: t_dep.append(t)
        elif data.EV_profile(t) == 0 and previous == True:
            t_dep.append(t)
            previous = False
        elif previous == False and data.EV_profile(t) == 0:
            previous = False
        else:
            print(t, "EV PROFILE EXTRACTION ERROR")             
    
    return EV_initial_SOC, t_arr, t_dep 

def configure_solver(model, solver_name="gurobi"):
    """
    :param model: the model to be solved
    :param solver_name: solver name as a string
    :return: a configured pyomo solver object
    """
    opt = SolverFactory('gurobi_persistent')
    opt.set_instance(model)
    return opt


def export_model(model):
    """
    :param model: Pyomo optimization model
    :return: Nothing, dumps lp file sizing_opt.lp
    """
    model.write(filename="sizing_opt.lp",
                format=ProblemFormat.cpxlp,
                io_options={"symbolic_solver_labels": True})
 
import plotly.graph_objects as go   
def plot_and_print(model, data, sizing=False, plot_soc=True):
    delta_t=RESOLUTION_IN_MINUTES/60
    conso = np.array([-data.consumption(p) for p in model.periods])
    PV =  np.array([model.gen_PV_sp[p].value for p in model.periods])
    steer =  np.array([model.steer_genset_sp[p].value for p in model.periods])
    imp =  np.array([model.imp_grid[p].value for p in model.periods])
    exp =  np.array([-model.exp_grid[p].value for p in model.periods])
    charge =  np.array([-model.charge_storage_sp[p].value for p in model.periods])
    discharge =  np.array([model.discharge_storage_sp[p].value for p in model.periods])
    if model.storage_energy_capacity.value > 0:
        soc =  np.array([model.soc[p].value/(model.storage_energy_capacity.value) for p in model.periods])
    else:
        soc = np.array([0 for p in model.periods])
        
    if BONUS:
        EV_charge =  np.array([-model.charge_EV_sp[p].value for p in model.periods])
        EV_discharge =  np.array([model.discharge_EV_sp[p].value for p in model.periods])
        EV_soc =  np.array([model.EV_soc[p].value/EV_CAPACITY for p in model.periods])
        
    fig, ax1 = plt.subplots(1, 1, dpi=720)
    if BONUS:
        1
        #ax1.stackplot(model.periods, conso,exp, charge, EV_charge, labels=['Consumption','Grid export', 'Battery charge', 'EV charge'])
        #ax1.stackplot(model.periods,PV,steer,imp,discharge,EV_discharge, labels=['PV','Genset','Grid import', 'Battery discharge', 'EV discharge'])
    else:
        ax1.stackplot(model.periods, conso,exp, charge, labels=['Consumption','Grid export', 'Battery charge'])
        #ax1.stackplot(model.periods ,charge, labels=[ 'Battery charge'])
        ax1.stackplot(model.periods,PV,steer,imp,discharge, labels=['PV','Genset','Grid import', 'Battery discharge'])
        #ax1.stackplot(model.periods,discharge, labels=['Battery discharge'])
    if plot_soc:
        ax2 = ax1.twinx()
        if BONUS:
            ax2.plot(model.periods, EV_soc)#, '.', markersize=2)
        ax2.plot(model.periods, soc)#, '.', markersize=2)
        ax2.set_ylabel('State of charge')
    

    ax1.legend(fontsize='small',ncol=2)

    if sizing:
        print("\nSizing decisions:")
        print("  - Storage capacity: %d kWh"  % (
            model.storage_energy_capacity.value/1e3))
        print("  - Storage inverter capacity: %d kVA"  % (
            model.storage_inverter_capacity.value/1e3))
        print("  - PV capacity: %d kWp"  % (model.PV_capacity.value/1e3))
        print("  - PV inverter capacity: %d kVA"  % (
            model.PV_inverter_capacity.value/1e3))
        if GRID_CONNECTED: print("  - Grid connection capacity: %d A"  % (model.grid_capacity.value))
        else: print("Grid connection capacity: 0 A")
        print("  - Genset capacity: %d kVA"  % (model.genset_capacity.value/1e3))
        
    
    print("\nTotal negative energy (consumption):")
    print("  - Base consumption: %d kWh" % (sum(conso)*delta_t/1e3))
    if GRID_CONNECTED: print("  - Grid export: %d kWh" % (sum(exp)*delta_t/1e3))
    print("  - Battery charge: %d kWh" % (sum(charge)*delta_t/1e3))
    if BONUS: print("  - EV charge: %d kWh" % (sum(EV_charge)*delta_t/1e3))
    
    
    print("\nTotal positive energy (production):")
    print("  - PV production: %d kWh" % (sum(PV)*delta_t/1e3))
    print("  - Genset production: %d kWh" % (sum(steer)*delta_t/1e3))
    if GRID_CONNECTED: print("  - Grid import: %d kWh" % (sum(imp)*delta_t/1e3))
    print("  - Battery discharge: %d kWh"% (sum(discharge)*delta_t/1e3))
    if BONUS: print("  - EV discharge: %d kWh"% (sum(EV_discharge)*delta_t/1e3))
    
    print("\nCost for the considered period [€]: %d" %(value(model.obj)))
    
    #print("  - Self Sufficiency: %d percentage"  % ((sum(conso)*delta_t/1e3)-(sum(imp)*delta_t/1e3))/(sum(conso)*delta_t/1e3))
    print("Performance Metrics")
    self_sufficiency_percentage = 100 * ( (sum((conso) * delta_t / 1e3)+ (sum((imp) * delta_t / 1e3))) / (sum((conso) * delta_t / 1e3)))
    print("   -Self Sufficiency: %.2f Percentage" % self_sufficiency_percentage)

    self_consumption_percentage = 100 * ( (sum((conso) * delta_t / 1e3)- (sum((exp) * delta_t / 1e3))) / (sum((conso) * delta_t / 1e3)))
    print("   -Self Consumption: %.2f Percentage" % self_consumption_percentage)
    
    variable_emissions = FUEL_CO2*(sum((steer) * delta_t / 1e3)) + GRID_CO2*(sum((imp) * delta_t / 1e3))-GRID_CO2*(sum((exp) * delta_t / 1e3))
    fixed_emissions =  (PV_CO2*model.PV_capacity.value + STORAGE_CO2*model.storage_energy_capacity.value + GENSET_CO2*model.genset_capacity.value + INVERTER_CO2*(model.storage_inverter_capacity.value+model.PV_inverter_capacity.value))/INVESTMENT_HORIZON
    total_emissions = variable_emissions+fixed_emissions
    
    print("   -Variable Emsissions %.2f grams" %  variable_emissions)
    print("   -Fixed Emissions: %.2f grams" %  fixed_emissions)
    print("   -Total Yearly Emissions: %.2f grams" %  total_emissions)

    return steer, PV, imp, exp , charge, discharge