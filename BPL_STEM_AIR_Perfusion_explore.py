# setup application functions BPL_STEM_AIR_Perfuson, dependent on previous import of functions from fmu_explore 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-28 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

def newplot(title='Stem cell perfusion cultivation', plotType='TimeSeries'):
   """ Standard plot window,
        title = '' """
       
   # Reset pens
   resetPen()
       
   # Plot diagram 
   if plotType == 'Basic':
      
      ax1 = plt.subplot(6,1,1)
      ax2 = plt.subplot(6,1,2)
      ax3 = plt.subplot(6,1,3)
      ax4 = plt.subplot(6,1,4)
      ax5 = plt.subplot(6,1,5)
      ax6 = plt.subplot(6,1,6)
      
      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
      ax.append(ax4)
      ax.append(ax5)
      ax.append(ax6)

      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('N [1E6]')
            
      ax[1].grid()
      ax[1].set_ylabel('G [g/L]')

      ax[2].grid()
      ax[2].set_ylabel('L [g/L]')
      
      ax[3].grid()
      ax[3].set_ylabel('DO [mg/L]')
      
      ax[4].grid()
      ax[4].set_ylabel('Vcc [L]')

      ax[5].grid()
      ax[5].set_ylabel('F [L/h]')
      ax[5].set_xlabel('Time [h]')
        
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['N'], color='b', linestyle=linetype)")       
      diagrams.append("ax[1].plot(t,sim_res['G'], color='b', linestyle=linetype)")  
      diagrams.append("ax[1].set_ylim(0)") 
      diagrams.append("ax[2].plot(t,sim_res['L'], color='b', linestyle=linetype)")  
      diagrams.append("ax[3].plot(t,sim_res['DO'], color='b', linestyle=linetype)")  
      diagrams.append("ax[3].set_ylim(0)") 
      diagrams.append("ax[4].plot(t,sim_res['Vcc'], color='b', linestyle=linetype)")  
      diagrams.append("ax[4].set_ylim([0, 0.050])")    
      diagrams.append("ax[5].step(t,sim_res['F'], color='b', linestyle=linetype)") 

   elif plotType == 'Comprehensive':

      ax11 = plt.subplot(5,2,1); ax12 = plt.subplot(5,2,2)
      ax21 = plt.subplot(5,2,3); ax22 = plt.subplot(5,2,4)
      ax31 = plt.subplot(5,2,5); ax32 = plt.subplot(5,2,6)
      ax41 = plt.subplot(5,2,7); ax42 = plt.subplot(5,2,8)
      ax51 = plt.subplot(5,2,9); ax52 = plt.subplot(5,2,10)
      
      ax.clear()
      ax.append(ax11) # 0
      ax.append(ax12) # 1
      ax.append(ax21) # 2
      ax.append(ax22) # 3
      ax.append(ax31) # 4
      ax.append(ax32) # 5
      ax.append(ax41) # 6
      ax.append(ax42) # 7
      ax.append(ax51) # 8
      ax.append(ax52) # 9
      
      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('N [1E6]')

      ax[2].grid()
      ax[2].set_ylabel('G [g/L]')

      ax[4].grid()
      ax[4].set_ylabel('L [g/L]')
      
      ax[6].grid()
      ax[6].set_ylabel('DO [mg/L]')

      ax[8].grid()
      ax[8].set_ylabel('F [L/h]')
      ax[8].set_xlabel('Time [h]')

      ax[1].grid()
      ax[1].set_ylabel('qN [1/h]')
      
      ax[3].grid()
      ax[3].set_ylabel('qG [g/(h*1E6)]')

      ax[5].grid()
      ax[5].set_ylabel('qL [g/(h*1E6)]')

      ax[7].grid()
      ax[7].set_ylabel('qO2 [mg/(h*1E6)]')

      ax[9].grid()
      ax[9].set_ylabel('OUR [mg/h]')
      ax[9].set_xlabel('Time [h]')


      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(t,sim_res['N'], color='b', linestyle=linetype)")       
      diagrams.append("ax[2].plot(t,sim_res['G'], color='b', linestyle=linetype)")  
      diagrams.append("ax[2].set_ylim(0)") 
      diagrams.append("ax[4].plot(t,sim_res['L'], color='b', linestyle=linetype)")  
      diagrams.append("ax[6].plot(t,sim_res['DO'], color='b', linestyle=linetype)")  
      diagrams.append("ax[6].set_ylim(0)")     
      diagrams.append("ax[8].step(t,sim_res['F'], color='b', linestyle=linetype)")       

      diagrams.append("ax[1].plot(t,sim_res['qN'], color='b', linestyle=linetype)") 
      diagrams.append("ax[1].step(t,sim_res['qNmax'], color='r', linestyle=linetype)") 
      diagrams.append("ax[3].plot(t,sim_res['qG'], color='b', linestyle=linetype)") 
      diagrams.append("ax[5].plot(t,sim_res['qLp']-sim_res['qLc'], color='b', linestyle=linetype)") 
      diagrams.append("ax[5].plot(t,sim_res['qLp'], color='g', linestyle=linetype)") 
      diagrams.append("ax[5].step(t,sim_res['qLc'], color='r', linestyle=linetype)") 
      diagrams.append("ax[7].step(t,sim_res['qO2'], color='r', linestyle=linetype)") 
      diagrams.append("ax[9].plot(t,sim_res['N']*sim_res['qO2'], color='b', linestyle=linetype)") 
             
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
#  Startup
#------------------------------------------------------------------------------------------------------------------

FMU_explore_info()