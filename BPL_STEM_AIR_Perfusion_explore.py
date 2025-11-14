# Figure - Simulation of perfusion reactor with STEM culture
#          with functions added to facilitate explorative simulation work 
#
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2024-08-31 - Created
# 2024-09-11 - Changed name of FMU
# 2024-10-03 - Changed name of FMU
# 2025-07-22 - Updated to MSL 4.1.0 and just the information text
# 2025-11-10 - Update FMU-explore 1.0.2 
# 2025-11-13 - Removed global declaration outside the funtions
# 2025-11-14 - FMU-explore 1.0.2 corrected
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.image as img
import zipfile 

from pyfmi import load_fmu
from pyfmi.fmi import FMUException

from itertools import cycle
from importlib_metadata import version   # included in future Python 3.8

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

# Simulation time
simulationTime = 1000.0
prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

# Define process_diagram
fmu_process_diagram = 'BPL_GUI_STEM_Perfusion_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------
   
# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = model.get_states_list()
stateValue.update(timeDiscreteStates)

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

def newplot(title='Stem cell perfusion cultivation', plotType='TimeSeries'):
   """ Standard plot window,
        title = '' """
    
   # Transfer of argument to global variable
   global ax1, ax2, ax3, ax4, ax5, ax6
   global ax11, ax21, ax31, ax41, ax51
   global ax12, ax22, ax32, ax42, ax52
   
   # Reset pens
   setLines()
       
   # Plot diagram 
   if plotType == 'Basic':
      
      plt.figure()
      ax1 = plt.subplot(6,1,1)
      ax2 = plt.subplot(6,1,2)
      ax3 = plt.subplot(6,1,3)
      ax4 = plt.subplot(6,1,4)
      ax5 = plt.subplot(6,1,5)
      ax6 = plt.subplot(6,1,6)

      ax1.set_title(title)
      ax1.grid()
      ax1.set_ylabel('N [1E6]')
            
      ax2.grid()
      ax2.set_ylabel('G [g/L]')

      ax3.grid()
      ax3.set_ylabel('L [g/L]')
      
      ax4.grid()
      ax4.set_ylabel('DO [mg/L]')
      
      ax5.grid()
      ax5.set_ylabel('Vcc [L]')

      ax6.grid()
      ax6.set_ylabel('F [L/h]')
      ax6.set_xlabel('Time [h]')
        
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(t,sim_res['N'], color='b', linestyle=linetype)")       
      diagrams.append("ax2.plot(t,sim_res['G'], color='b', linestyle=linetype)")  
      diagrams.append("ax2.set_ylim(0)") 
      diagrams.append("ax3.plot(t,sim_res['L'], color='b', linestyle=linetype)")  
      diagrams.append("ax4.plot(t,sim_res['DO'], color='b', linestyle=linetype)")  
      diagrams.append("ax4.set_ylim(0)") 
      diagrams.append("ax5.plot(t,sim_res['Vcc'], color='b', linestyle=linetype)")  
      diagrams.append("ax5.set_ylim([0, 0.050])")    
      diagrams.append("ax6.step(t,sim_res['F'], color='b', linestyle=linetype)") 

   elif plotType == 'Comprehensive':

      plt.figure()
      ax11 = plt.subplot(5,2,1); ax12 = plt.subplot(5,2,2)
      ax21 = plt.subplot(5,2,3); ax22 = plt.subplot(5,2,4)
      ax31 = plt.subplot(5,2,5); ax32 = plt.subplot(5,2,6)
      ax41 = plt.subplot(5,2,7); ax42 = plt.subplot(5,2,8)
      ax51 = plt.subplot(5,2,9); ax52 = plt.subplot(5,2,10)

      ax11.set_title(title)
      ax11.grid()
      ax11.set_ylabel('N [1E6]')

      ax21.grid()
      ax21.set_ylabel('G [g/L]')

      ax31.grid()
      ax31.set_ylabel('L [g/L]')
      
      ax41.grid()
      ax41.set_ylabel('DO [mg/L]')

      ax51.grid()
      ax51.set_ylabel('F [L/h]')
      ax51.set_xlabel('Time [h]')

      ax12.grid()
      ax12.set_ylabel('qN [1/h]')
      
      ax22.grid()
      ax22.set_ylabel('qG [g/(h*1E6)]')

      ax32.grid()
      ax32.set_ylabel('qL [g/(h*1E6)]')

      ax42.grid()
      ax42.set_ylabel('qO2 [mg/(h*1E6)]')

      ax52.grid()
      ax52.set_ylabel('OUR [mg/h]')
      ax52.set_xlabel('Time [h]')


      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax11.plot(t,sim_res['N'], color='b', linestyle=linetype)")       
      diagrams.append("ax21.plot(t,sim_res['G'], color='b', linestyle=linetype)")  
      diagrams.append("ax21.set_ylim(0)") 
      diagrams.append("ax31.plot(t,sim_res['L'], color='b', linestyle=linetype)")  
      diagrams.append("ax41.plot(t,sim_res['DO'], color='b', linestyle=linetype)")  
      diagrams.append("ax41.set_ylim(0)")     
      diagrams.append("ax51.step(t,sim_res['F'], color='b', linestyle=linetype)")       

      diagrams.append("ax12.plot(t,sim_res['qN'], color='b', linestyle=linetype)") 
      diagrams.append("ax12.step(t,sim_res['qNmax'], color='r', linestyle=linetype)") 
      diagrams.append("ax22.plot(t,sim_res['qG'], color='b', linestyle=linetype)") 
      diagrams.append("ax32.plot(t,sim_res['qLp']-sim_res['qLc'], color='b', linestyle=linetype)") 
      diagrams.append("ax32.plot(t,sim_res['qLp'], color='g', linestyle=linetype)") 
      diagrams.append("ax32.step(t,sim_res['qLc'], color='r', linestyle=linetype)") 
      diagrams.append("ax42.step(t,sim_res['qO2'], color='r', linestyle=linetype)") 
      diagrams.append("ax52.plot(t,sim_res['N']*sim_res['qO2'], color='b', linestyle=linetype)") 
             
def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""

   if name == 'culture':
      print('Reactor culture human-induced pluripotent stem cells - hiPSCs') 

   elif name in ['broth', 'liquidphase', 'liquid-phase''media']:

      X = model.get('liquidphase.X')[0]; X_description = model.get_variable_description('liquidphase.X'); X_mw = model.get('liquidphase.mw[1]')[0]
      G = model.get('liquidphase.G')[0]; G_description = model.get_variable_description('liquidphase.G'); G_mw = model.get('liquidphase.mw[2]')[0]
      L = model.get('liquidphase.L')[0]; L_description = model.get_variable_description('liquidphase.L'); L_mw = model.get('liquidphase.mw[3]')[0]

      print('Reactor broth substances included in the model')
      print()
      print(X_description, '   index = ', X, 'molecular weight = ', X_mw, 'Da')
      print(G_description, ' index = ', G, 'molecular weight = ', G_mw, 'Da')
      print(L_description, ' index = ', L, 'molecular weight = ', L_mw, 'Da')

   elif name in ['parts']:
      describe_parts(component_list_minimum)

   elif name in ['MSL']:
      describe_MSL()

   else:
      describe_general(name, decimals)

#------------------------------------------------------------------------------------------------------------------
#  General code 
FMU_explore = 'FMU-explore version 1.0.2'
#------------------------------------------------------------------------------------------------------------------

# Define function par() for parameter update
def par(*x, parValue=parValue, **x_kwarg):
   """ Set parameter values if available in the predefined dictionaryt parValue. """
   x_kwarg.update(*x)
   x_temp = {}
   for key in x_kwarg.keys():
      if key in parValue.keys():
         x_temp.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an accessible parameter - check the spelling')
   parValue.update(x_temp)
   
   parErrors = [requirement for requirement in parCheck if not(eval(requirement))]
   if not parErrors == []:
      print('Error - the following requirements do not hold:')
      for index, item in enumerate(parErrors): print(item)

# Define function init() for initial values update
def init(*x, parValue=parValue, **x_kwarg):
   """ Set initial values and the name should contain string '_start' to be accepted.
       The function can handle general parameter string location names if entered as a dictionary. """
   x_kwarg.update(*x)
   x_init={}
   for key in x_kwarg.keys():
      if '_start' in key: 
         x_init.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an initial value, use par() instead - check the spelling')
   parValue.update(x_init)

# Define how to read dictionary for parameter values
def readParValue(file, sheet, parValue=parValue):
   """ Read parameter short names and values from an Excel-file from defined sheet. For use in the notebook!
       Return a dictionary."""
   parValue_local = {} 
   table = pd.ExcelFile(file).parse(sheet)
   for k in list(range(len(table))):
      parValue_local[table['Par'][k]] = table['Value'][k]
   parValue.update(parValue_local)

# Define how to read dictionary for parameter location
def readParLocation(file, parLocation=parLocation):
   """ Read parameter short and long names from an Excel-file sheet by sheet. For use in the notebook!
       Return a dictionary."""
   sheets = ['initial_values','feed_AB', 'feed_G', 'culture', 'broth_decay']
   parLocation_local = {}
   for sheet in sheets:
      table = pd.ExcelFile(file).parse(sheet)
      for k in list(range(len(table))):
         parLocation_local[table['Par'][k]] = table['Location'][k]
   parLocation.update(parLocation_local)
      
def disp(name='', decimals=3, mode='short', parValue=parValue, parLocation=parLocation):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   global model

   def dict_reverser(d):
      seen = set()
      return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parValue.keys()]:
         if name in Location:
            if type(model.get(Location)[0]) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model.get(Location)[0],decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model.get(Location)[0])               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parValue.keys():
            if name in parName:
               if type(model.get(Location)[0]) != np.bool_:
                  print(parName,':', np.round(model.get(parLocation[parName])[0],decimals))
               else: 
                  print(parName,':', model.get(parLocation[parName])[0])
   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parValue.keys()]:
         if name in Location:
            if type(model.get(Location)[0]) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model.get(Location)[0],decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parValue.keys():
            if name in parName:
               if type(model.get(Location)[0]) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model.get(parLocation[parName])[0],decimals))

# Line types
def setLines(lines=['-','--',':','-.']):
   """Set list of linetypes used in plots"""
   global linecycler
   linecycler = cycle(lines)

# Show plots from sim_res, just that
def show(diagrams=diagrams):
   """Show diagrams chosen by newplot()"""
   # Plot pen
   linetype = next(linecycler)    
   # Plot diagrams 
   for command in diagrams: eval(command)

# Simulation
def simu(simulationTimeLocal=simulationTime, mode='Initial', options=opts_std, \
         diagrams=diagrams,timeDiscreteStates=timeDiscreteStates, stateValue=stateValue, \
         parValue=parValue, parLocation=parLocation, fmu_model=fmu_model):         
   """Model loaded and given intial values and parameter before,
      and plot window also setup before."""
    
   # Global variables
   global model, prevFinalTime, sim_res, t
   
   # Simulation flag
   simulationDone = False
   
   # Transfer of argument to global variable
   simulationTime = simulationTimeLocal 
      
   # Check parValue
   value_missing = 0
   for key in parValue.keys():
      if parValue[key] in [np.nan, None, '']:
         print('Value missing:', key)
         value_missing =+1
   if value_missing>0: return
         
   # Load model
   if model is None:
      model = load_fmu(fmu_model) 
   model.reset()
      
   # Run simulation
   if mode in ['Initial', 'initial', 'init']:
      # Set parameters and intial state values:
      for key in parValue.keys():
         model.set(parLocation[key],parValue[key])   
      # Simulate
      sim_res = model.simulate(final_time=simulationTime, options=options)  
      simulationDone = True
   elif mode in ['Continued', 'continued', 'cont']:

      if prevFinalTime == 0: 
         print("Error: Simulation is first done with default mode = init'")      
      else:
         
         # Set parameters and intial state values:
         for key in parValue.keys():
            model.set(parLocation[key],parValue[key])                

         for key in stateValue.keys():
            if not key[-1] == ']':
               if key[-3:] == 'I.y': 
                  model.set(key[:-10]+'I_start', stateValue[key]) 
               elif key[-3:] == 'D.x': 
                  model.set(key[:-10]+'D_start', stateValue[key]) 
               else:
                  model.set(key+'_start', stateValue[key])
            elif key[-3] == '[':
               model.set(key[:-3]+'_start'+key[-3:], stateValue[key]) 
            elif key[-4] == '[':
               model.set(key[:-4]+'_start'+key[-4:], stateValue[key]) 
            elif key[-5] == '[':
               model.set(key[:-5]+'_start'+key[-5:], stateValue[key]) 
            else:
               print('The state vecotr has more than 1000 states')
               break

         # Simulate
         sim_res = model.simulate(start_time=prevFinalTime,
                                 final_time=prevFinalTime + simulationTime,
                                 options=options) 
         simulationDone = True             
   else:
      print("Simulation mode not correct")

   if simulationDone:
    
      # Extract data
      t = sim_res['time']
 
      # Plot diagrams
      linetype = next(linecycler)    
      for command in diagrams: eval(command)
            
      # Store final state values stateValue:
      for key in list(stateValue.keys()): stateValue[key] = model.get(key)[0]        

      # Store time from where simulation will start next time
      prevFinalTime = model.time
   
   else:
      print('Error: No simulation done')
      
# Describe model parts of the combined system
def describe_parts(component_list=[]):
   """List all parts of the model""" 
       
   def model_component(variable_name):
      i = 0
      name = ''
      finished = False
      if not variable_name[0] == '_':
         while not finished:
            name = name + variable_name[i]
            if i == len(variable_name)-1:
                finished = True 
            elif variable_name[i+1] in ['.', '(']: 
                finished = True
            else: 
                i=i+1
      if name in ['der', 'temp_1', 'temp_2', 'temp_3', 'temp_4', 'temp_5', 'temp_6', 'temp_7']: name = ''
      return name
    
   variables = list(model.get_model_variables().keys())
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))
   
def describe_MSL(flag_vendor=flag_vendor):
   """List MSL version and components used"""
   print('MSL:', MSL_usage)
 
# Describe parameters and variables in the Modelica code
def describe_general(name, decimals, parLocation=parLocation):
  
   if name == 'time':
      description = 'Time'
      unit = 'h'
      print(description,'[',unit,']')
      
   elif name in parLocation.keys():
      description = model.get_variable_description(parLocation[name])
      value = model.get(parLocation[name])[0]
      try:
         unit = model.get_variable_unit(parLocation[name])
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)            
      else:
        print(description, ':', np.round(value, decimals), '[',unit,']')
                  
   else:
      description = model.get_variable_description(name)
      value = model.get(name)[0]
      try:
         unit = model.get_variable_unit(name)
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)     
      else:
         print(description, ':', np.round(value, decimals), '[',unit,']')
         
# Plot process diagram
def process_diagram(fmu_model=fmu_model, fmu_process_diagram=fmu_process_diagram):   
   try:
       process_diagram = zipfile.ZipFile(fmu_model, 'r').open('documentation/processDiagram.png')
   except KeyError:
       print('No processDiagram.png file in the FMU, but try the file on disk.')
       process_diagram = fmu_process_diagram
   try:
       plt.imshow(img.imread(process_diagram))
       plt.axis('off')
       plt.show()
   except FileNotFoundError:
       print('And no such file on disk either')
         
# Describe framework
def BPL_info():
   print()
   print('Model for the process has been setup. Key commands:')
   print(' - par()       - change of parameters and initial values')
   print(' - init()      - change initial values only')
   print(' - simu()      - simulate and plot')
   print(' - newplot()   - make a new plot')
   print(' - show()      - show plot from previous simulation')
   print(' - disp()      - display parameters and initial values from the last simulation')
   print(' - describe()  - describe culture, broth, parameters, variables with values/units')
   print()
   print('Note that both disp() and describe() takes values from the last simulation')
   print('and the command process_diagram() brings up the main configuration')
   print()
   print('Brief information about a command by help(), eg help(simu)') 
   print('Key system information is listed with the command system_info()')

def system_info():
   """Print system information"""
   FMU_type = model.__class__.__name__
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -PyFMI:', version('pyfmi'))
   print(' -FMU by:', model.get_generation_tool())
   print(' -FMI:', model.get_version())
   print(' -Type:', FMU_type)
   print(' -Name:', model.get_name())
   print(' -Generated:', model.get_generation_date_and_time())
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)
   
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()