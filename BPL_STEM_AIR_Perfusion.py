# setup data STEM_AIR_Perfusion
# Author: Jan Peter Axelsson

#------------------------------------------------------------------------------------------------------------------
# 2026-08-28 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt 
from pyfmi import load_fmu

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   flag_vendor = 'JM'
   flag_type = 'CS'
   fmu_model ='BPL_STEM_AIR_Reactor_windows_jm_cs.fmu'        
   model = load_fmu(fmu_model, log_level=0)  
elif platform.system() == 'Linux':  
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-compiled OpenModelica') 
      if flag_type in ['CS','cs']:         
         fmu_model =''    
         model = load_fmu(fmu_model, log_level=0) 
      if flag_type in ['ME','me']:         
         fmu_model = 'BPL_STEM_AIR_Reactor_linux_om_me.fmu'     
         model = load_fmu(fmu_model, log_level=0)
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = model.simulate_options()
   opts_std['silent_mode'] = True
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary'     
elif flag_type in ['ME', 'me']:
   opts_std = model.simulate_options()
   opts_std["CVode_options"]["verbosity"] = 50 
   opts_std['ncp'] = 500 
   opts_std['result_handling'] = 'binary'  
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm', 'OM', 'om']:
   MSL_usage = '4.1.0 - used components: RealInput, RealOutput, CombiTimeTable, Types' 
   MSL_version = '4.1.0'
   BPL_version = 'BPL - not used' 
   
#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, ax, lines
#------------------------------------------------------------------------------------------------------------------

# Simulation time
simulationTime = 1000.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = model.get_states_list()
stateValue.update(timeDiscreteStates)

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

# Define process_diagram
fmu_process_diagram = 'BPL_GUI_STEM_Perfusion_process_diagram_om.png'

# Create parValue
parValue = {}
parValue['Vcc'] = 0.040                # L
parValue['N_start'] = 50               # 1E6 
parValue['DO_start'] = 100             # mg/L
parValue['qm'] = 1.0e-6                # mg/(h*1E6)
parValue['Yns'] = 0.013                # mg/1E6
parValue['qLpmax'] = 0.008             #
parValue['scale'] = 1000.0             # Correction of L equation Vcc
parValue['CL0'] = 8.0                  # mg/L
parValue['OTR'] = 21                   # mg/h  

parLocation = {}
parLocation['Vcc'] = 'Vcc'
parLocation['N_start'] = 'N_start'
parLocation['DO_start'] = 'DO_start'
parLocation['qm'] = 'qm'
parLocation['Yns'] = 'Yns'
parLocation['qLpmax'] = 'qLpmax'
parLocation['scale'] = 'scale'
parLocation['CL0'] = 'CL0'
parLocation['OTR'] = 'OTR'

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []

# Create list of diagrams to be plotted by simu()
diagrams = []

# Create an empty list axes to be defined in newplot() and plotted by simu() or show()
ax = []

# Create list of pens for the diagrams
lines = ['-','--',':','-.']