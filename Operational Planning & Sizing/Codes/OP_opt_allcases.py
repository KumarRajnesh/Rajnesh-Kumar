from pyomo.core.base import Var, Constraint, Objective, NonNegativeReals, Reals, Binary, Integers, maximize, Param
import pyomo.environ as pyo
from pyomo.core.kernel import value
from pyomo.opt import SolverStatus, TerminationCondition
from pyomo.core.base import ConcreteModel, RangeSet, Set
import utils
from utils import extract_EV_data

from parameters import N_PERIODS, START_DATE, RESOLUTION_IN_MINUTES, N_DAYS, INVESTMENT_HORIZON
from parameters import BONUS
from parameters import GRID_VOLTAGE, ANNUAL_CONSUMPTION
from parameters import CHARGE_EFFICIENCY, DISCHARGE_EFFICIENCY, STORAGE_MAX_P_RATE,PROSUMER_TARIF, PV_CAPACITY_PRICE
from parameters import FUEL_PRICE_COEFF, GRID_EXPORT_PRICE, GRID_IMPORT_PRICE, GRID_CAPACITY_PRICE, INITIAL_SOC

# Load data                                             
data = utils.SizingData(start_date=START_DATE, yearly_cons=ANNUAL_CONSUMPTION)

# Create the pyomo model
model = ConcreteModel()

# Create sets
# Time periods of length RESOLUTION_IN_MINUTES
model.periods = RangeSet(N_PERIODS)                                             # Goes from 1 to N_PERIODS

# Create variables
# Operation variables
model.soc = Var(model.periods, within=NonNegativeReals)                         # Storage state of charge [Wh]
model.charge_storage_sp = Var(model.periods, within=NonNegativeReals)           # Storage charge setpoint [W]
model.discharge_storage_sp = Var(model.periods, within=NonNegativeReals)        # Storage discharge setpoint [W]

model.EV_soc = Var(model.periods, within=NonNegativeReals)                      # EV state of charge [Wh]
model.charge_EV_sp = Var(model.periods, within=NonNegativeReals)                # EV charge setpoint [W]
model.discharge_EV_sp = Var(model.periods, within=NonNegativeReals)             # EV discharge setpoint [W]

model.gen_PV_sp = Var(model.periods, within=NonNegativeReals)                   # PV generation [W]
model.imp_grid = Var(model.periods, within=NonNegativeReals)                    # [W]
model.exp_grid = Var(model.periods, within=NonNegativeReals)                    # [W]
model.steer_genset_sp = Var(model.periods, within=NonNegativeReals)             # genset setpoint [W]

# Sizing PARAMETERS
model.storage_energy_capacity = Param(default=45e3)                             # Battery capacity is 45 kWh
model.storage_inverter_capacity = Param(default=10e3)                           # W
model.PV_capacity = Param(default=10e3)                                         # PV system capacity is 10 kWp
model.PV_inverter_capacity = Param(default=10e3)                                # Pnom of the inverter is 10 kW
model.grid_capacity = Param(default=32)                                         # Connecion limit in [A]
model.genset_capacity = Param(default=10e3)                                     # Genset Pnom is 10 kW


# Create constraints:                                   
#           use data.consumption(p) to get consumption for timestep p
#               data.PV_generation(p) to get Pmax [pu] for timestep p
#               data.EV_connected(p) to get connection status of the EV (0 or 1)
#           p should go from 1 to N_PERIODS

def power_balance_rule(m, p):
    
    prod_expression = m.gen_PV_sp[p] + m.steer_genset_sp[p] +  m.discharge_storage_sp[p] + m.imp_grid[p] + m.discharge_EV_sp[p]
    cons_expression = data.consumption(p) + m.exp_grid[p] + m.charge_EV_sp[p] + m.charge_storage_sp[p]

    return prod_expression == cons_expression                                   

def storage_level_rule(m, p):
    if p==1:
        return m.soc[p]==INITIAL_SOC*m.storage_energy_capacity
    else:
        return m.soc[p] == m.soc[p-1] + (m.charge_storage_sp[p] * CHARGE_EFFICIENCY - m.discharge_storage_sp[p] / DISCHARGE_EFFICIENCY) * RESOLUTION_IN_MINUTES / 60
        
def storage_max_capacity_rule(m, p):
    return m.soc[p] <= m.storage_energy_capacity*0.9
def storage_min_capacity_rule(m, p):
    return m.soc[p] >= m.storage_energy_capacity*0.2  #minimum SOC requirement
    #TODO
def storage_charge_power_rule1(m, p):
    return m.charge_storage_sp[p] <= STORAGE_MAX_P_RATE*m.storage_energy_capacity

def storage_discharge_power_rule1(m, p):
    return m.discharge_storage_sp[p] <= STORAGE_MAX_P_RATE*m.storage_energy_capacity 

def storage_charge_power_rule2(m, p):
    return m.charge_storage_sp[p] <= m.storage_inverter_capacity 

def storage_discharge_power_rule2(m, p):
    return m.discharge_storage_sp[p] <= m.storage_inverter_capacity 

def PV_generation_rule1(m, p):
    return m.gen_PV_sp[p] <= m.PV_inverter_capacity

def PV_generation_rule2(m, p):
    return m.gen_PV_sp[p] <= m.PV_capacity*data.PV_generation(p)

def import_limit_rule1(m, p):
    return m.imp_grid[p] <= m.grid_capacity*GRID_VOLTAGE

def import_limit_rule2(m, p):
    return m.imp_grid[p] <= data.consumption(p) #TODO

def export_limit_rule(m, p):
    return m.exp_grid[p] <= m.grid_capacity*GRID_VOLTAGE

def genset_generation_rule(m, p):
    return m.steer_genset_sp[p] <= m.genset_capacity


model.power_balance_cstr = Constraint(model.periods, rule=power_balance_rule)

model.storage_level_cstr = Constraint(model.periods, rule=storage_level_rule)
model.storage_min_capacity_cstr = Constraint(model.periods, rule=storage_min_capacity_rule)
model.storage_max_capacity_cstr = Constraint(model.periods, rule=storage_max_capacity_rule)
model.storage_charge_power_cstr1 = Constraint(model.periods, rule=storage_charge_power_rule1)
model.storage_discharge_power_cstr1 = Constraint(model.periods, rule=storage_discharge_power_rule1)
model.storage_charge_power_cstr2 = Constraint(model.periods, rule=storage_charge_power_rule2)
model.storage_discharge_power_cstr2 = Constraint(model.periods, rule=storage_discharge_power_rule2)

model.PV_generation_cstr1 = Constraint(model.periods, rule=PV_generation_rule1)
model.PV_generation_cstr2 = Constraint(model.periods, rule=PV_generation_rule2)

model.import_limit_cstr1 = Constraint(model.periods, rule=import_limit_rule1)
model.import_limit_cstr2 = Constraint(model.periods, rule=import_limit_rule2)
model.export_limit_cstr = Constraint(model.periods, rule=export_limit_rule)

model.genset_generation_cstr = Constraint(model.periods, rule=genset_generation_rule)

# --------------------BONUS-------------------
if BONUS:
    from parameters import EV_TARGET_SOC, EV_CAPACITY, EV_INVERTER_CAPACITY, YEARLY_KM
    EV_initial_SOC, t_arr, t_dep = extract_EV_data(data, YEARLY_KM)             # Create t_init array, t_target and soc_init
    n_EV_con = len(t_arr)                                                       # Number of EV_connections
    
    # EV_initial_SOC: list of states of charge when the EV arrives [Wh] indices go from 0 to n_EV_con-1
    # t_arr: list of time indices for EV arrival. Numbers in this list can go from 1 to N_PERIODS, indices of this list go from 0 to n_EV_con-1
    #       you should ensure that at the begining of each timestep in this list, the EV_soc variable is equal to the corresponding EV_initial_SOC
    # t_dep: list of time indices for EV departure. Numbers in this list can go from 1 to N_PERIODS, indices of this list go from 0 to n_EV_con-1
    #       you should ensure that at the end of each timestep in this list, the EV_soc variable is equal to EV_TARGET_SOC
    model.connections = Set(initialize=range(n_EV_con))                         # Goes from 0 to (number of connections-1)

    
    model.EV_capacity = Param(default=EV_CAPACITY)                              # Capacity of EV battery is 60 kWh
    model.EV_inverter_capacity = Param(default=EV_INVERTER_CAPACITY)            # Pnom of the inverter is 10 kVA
    
    
    def EV_target_rule(m,c):
        return m.EV_soc[t_dep[c]] == EV_TARGET_SOC                      
    
    def EV_max_capacity_rule(m, p):
        return m.EV_soc[p] <=  0.9* m.EV_capacity                   
    
    def EV_min_capacity_rule(m, p):
        return m.EV_soc[p] >= 0.2*m.EV_capacity                   
    
    def EV_charge_power_rule(m, p):
        if data.EV_connected(p) == 1:
            return m.charge_EV_sp[p] <= m.EV_inverter_capacity 
        else:
            return m.charge_EV_sp[p]==0                   
    
    def EV_discharge_power_rule(m, p):
        if data.EV_connected(p) == 1:
            return m.discharge_EV_sp[p] <= m.EV_inverter_capacity
        else:
            return m.discharge_EV_sp[p]==0                                        
    
    def EV_level_rule(m, p):
        if p in t_arr:
           return m.EV_soc[p] == EV_initial_SOC[t_arr.index(p)]
        elif p in t_dep:

           return m.EV_soc[p-1] + m.charge_EV_sp[p] == m.EV_soc[p]
        elif p==1:
            return Constraint.Skip
        elif  data.EV_connected(p) == 0 and data.EV_connected(p)==0: 
            return Constraint.Skip
        else:
            return m.EV_soc[p] == m.EV_soc[p-1] + (m.charge_EV_sp[p]  - m.discharge_EV_sp[p])*RESOLUTION_IN_MINUTES/60 
    
    def EV_simulatneous_rule(m,p):
        #No simultaneous charge and discharge
        return m.charge_EV_sp[p]*m.discharge_EV_sp[p] == 0


    model.EV_level_cstr = Constraint(model.periods, rule=EV_level_rule)
    model.EV_min_capacity_cstr = Constraint(model.periods, rule=EV_min_capacity_rule)
    model.EV_max_capacity_cstr = Constraint(model.periods, rule=EV_max_capacity_rule)
    model.EV_charge_power_cstr = Constraint(model.periods, rule=EV_charge_power_rule)
    model.EV_discharge_power_cstr = Constraint(model.periods, rule=EV_discharge_power_rule)
    model.EV_target_cstr = Constraint(model.connections, rule=EV_target_rule)
# --------------------------------------------
else:
    def EV_charge_power_rule2(m, p):                                           # Null when no EV
        return m.charge_EV_sp[p] <= 0
    
    def EV_discharge_power_rule2(m, p):                                        # Null when no EV
        return m.discharge_EV_sp[p] <= 0
    model.EV_charge_power_cstr = Constraint(model.periods, rule=EV_charge_power_rule2)
    model.EV_discharge_power_cstr = Constraint(model.periods, rule=EV_discharge_power_rule2)
    

#Create objective

SG = True  # selling to grid
GA = False # giving it away
NM = False # net metering

if SG:

    def objective(m):
        
        opex = GRID_CAPACITY_PRICE*m.grid_capacity+sum(GRID_IMPORT_PRICE * m.imp_grid[p]/1000*RESOLUTION_IN_MINUTES / 60 for p in m.periods)-sum(GRID_EXPORT_PRICE * m.exp_grid[p]/1000*RESOLUTION_IN_MINUTES / 60 for p in m.periods)+sum(FUEL_PRICE_COEFF * m.steer_genset_sp[p]/1000*RESOLUTION_IN_MINUTES / 60 for p in m.periods)
        return opex
    
if GA:
    def objective(m):
        
        opex = GRID_CAPACITY_PRICE*m.grid_capacity+sum(GRID_IMPORT_PRICE * m.imp_grid[p]/1000*RESOLUTION_IN_MINUTES / 60 for p in m.periods)+sum(FUEL_PRICE_COEFF * m.steer_genset_sp[p]/1000*RESOLUTION_IN_MINUTES / 60 for p in m.periods)
        return opex


if NM:
    model.net_imports = Var(within=NonNegativeReals)
    
    def net_metering_rule(m,p):
            return m.net_imports == sum(m.imp_grid[p] for p in m.periods) - sum(m.exp_grid[p] for p in m.periods)
    
    model.net_metering_cstr = Constraint(rule=net_metering_rule)
    
    def objective(m):
        
        opex = PROSUMER_TARIF*m.PV_inverter_capacity+GRID_CAPACITY_PRICE*m.grid_capacity +sum(FUEL_PRICE_COEFF * m.steer_genset_sp[p]/1000*RESOLUTION_IN_MINUTES / 60 for p in m.periods) +GRID_IMPORT_PRICE*m.net_imports/1000*RESOLUTION_IN_MINUTES / 60
        return opex


model.obj = Objective(rule=objective, sense=pyo.minimize)


# Solve the model
print("Solving model...")
opt = utils.configure_solver(model)
results = opt.solve(model, tee=False, keepfiles=False)

# Check is the problem is feasible
status = results.solver.status
termination_condition = results.solver.termination_condition

print("\nSolver status:", status, termination_condition)
if status == 'ok' : utils.plot_and_print(model,data, plot_soc=False) # You can change plot_soc to False to improve readability
