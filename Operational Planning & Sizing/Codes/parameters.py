# This file contains parameters you may want to play with.

from datetime import datetime
# Parameters

## General simulation parameters
YEARLY_KM = 5000                                                                   # TODO: add your car usage [kM/year]
ANNUAL_CONSUMPTION = 7.5e6                                                         # TODO: Define your own annual energy consumption (base 4 MWh)
# Connection to the grid
GRID_CONNECTED = True                                                           # TODO: Change depending on the current case
BONUS = True                                                                   # TODO: Decide wether to model EV

N_DAYS = 365                                                              # Number of days to simulate, from start_date [day] #TODO: Change to 365 to finalize


# Time
START_DATE = datetime(2021, 1, 1, 0, 0, 0)                                      #TODO: You can change that to model other seasons, but there is only data for 2021
# Start date of the sizing run. No error handling if out of the data range.

RESOLUTION_IN_MINUTES = 15                                                      # Time step duration [min]
N_PERIODS = int(N_DAYS * 24 * 60 / RESOLUTION_IN_MINUTES)                       # Number of operation time steps
INVESTMENT_HORIZON = 20                                                         # [Years] Investment horizon

# Storage device sizing parameters
STORAGE_PRICE = 0.3                                                            # [EUR/Wh]
STORAGE_MAX_P_RATE = 1.8                                                        # [W/Wh]

# Storage device operational parameters
INITIAL_SOC = 0.2  # [/] * Battery capacity
CHARGE_EFFICIENCY = 0.95  # [/]
DISCHARGE_EFFICIENCY = 0.95  # [/]

# Inverter capacity
INVERTER_PRICE = 150e-3  # [EUR/VA]
PROSUMER_TARIF = 0.070 # [EUR/VA.year]

# PV
PV_CAPACITY_PRICE = 0.5                                                         # [EUR/Wp]
PV_MAX_CAPACITY = 10e3                                                           # TODO: Depending on your roofsize, you can relax afterwards and see effects


GRID_CAPACITY_PRICE = 3                                                        # Grid connection cost [EUR/A.year]
GRID_VOLTAGE = 230                                                              # [V]
GRID_IMPORT_PRICE = 0.4                                                         # [EUR/kWh]
GRID_EXPORT_PRICE = 0.04                                                   # [EUR/kWh]
# Genset
GENSET_CAPACITY_PRICE = 0.15                                                    # [€/VA]


FUEL_PRICE_COEFF = 2/(10.74*0.4)                                                # [EUR/kWh] cost/(energy_density*efficiency)

# EV
EV_CAPACITY = 60e3                                                              # [Wh]: Capacity of the EV battery
EV_TARGET_SOC = 0.8*EV_CAPACITY                                                 # [Wh]: SoC needed at disconnection of the EV.
EV_INVERTER_CAPACITY = 10e3                                                     # [W]: Nominal power of EV inverter

## CO2
# Variable emission
FUEL_CO2 = 0.300                                                                # [g/Wh]
GRID_CO2 = 0.174                                                                # [g/Wh]

# Fixed emission
PV_CO2 = 1.8e3 # g/Wp
STORAGE_CO2 = 150 #g/Wh
GENSET_CO2 = 10  #g/VA
INVERTER_CO2 = 60 #g/VA
