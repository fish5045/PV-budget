## Budgets 
# Version 1.0
# Author: Shang-En Li
# Last edit: 20220626


### 仍要處理的功能


#############
## imports

import numpy as np
import matplotlib.pyplot as plt
import json
import os
import datetime
from meteotools.grids import cylindrical_grid
from meteotools.fileprocess import get_wrfout_time
from multiprocessing import Pool


################
## functions

def plotting_PV_RZ(t):
    lhs = np.zeros([total_t, Nz, Nr])
    adv_a = np.zeros([total_t, Nz, Nr])
    adv_r = np.zeros([total_t, Nz, Nr])
    adv_z = np.zeros([total_t, Nz, Nr])
    dia_a = np.zeros([total_t, Nz, Nr])
    dia_r = np.zeros([total_t, Nz, Nr])
    dia_z = np.zeros([total_t, Nz, Nr])
    fric_a = np.zeros([total_t, Nz, Nr])
    fric_r = np.zeros([total_t, Nz, Nr])
    fric_z = np.zeros([total_t, Nz, Nr])
    for idx, i in enumerate(range(-total_t//2+1,total_t//2+1)):
        
        curr_time = get_wrfout_time(setting['start_yr'], setting['start_mo'], 
                                    setting['start_dy'], setting['start_hr'], 
                                    setting['start_mi'], setting['start_sc'],
                                    (i+t)*dt)
        cyl_curr = cylindrical_grid(f'./output/PV_budget_cyl_{curr_time}.nc')
        cyl_curr.ncfile_readvars('lhs', 'adv_r', 'adv_a', 'adv_z',
                                 'dia_r', 'dia_a', 'dia_z', 
                                 'fric_r', 'fric_a', 'fric_z')
        lhs[idx,:,:] = np.nanmean(cyl_curr.lhs, axis=1)
        adv_a[idx,:,:] = np.nanmean(cyl_curr.adv_a, axis=1)
        adv_r[idx,:,:] = np.nanmean(cyl_curr.adv_r, axis=1)
        adv_z[idx,:,:] = np.nanmean(cyl_curr.adv_z, axis=1)
        dia_a[idx,:,:] = np.nanmean(cyl_curr.dia_a, axis=1)
        dia_r[idx,:,:] = np.nanmean(cyl_curr.dia_r, axis=1)
        dia_z[idx,:,:] = np.nanmean(cyl_curr.dia_z, axis=1)
        fric_r[idx,:,:] = np.nanmean(cyl_curr.fric_r, axis=1)
        fric_a[idx,:,:] = np.nanmean(cyl_curr.fric_a, axis=1)
        fric_z[idx,:,:] = np.nanmean(cyl_curr.fric_z, axis=1)
        if i == 0:
            output_time = curr_time
            print(output_time)
    
    total_lhs = np.nanmean(lhs, axis=0)
    radial_adv = np.nanmean(adv_r, axis=0)
    azimuthal_adv = np.nanmean(adv_a, axis=0)
    vertical_adv = np.nanmean(adv_z, axis=0)
    diabat_r = np.nanmean(dia_r, axis=0)
    diabat_a = np.nanmean(dia_a, axis=0)
    diabat_z = np.nanmean(dia_z, axis=0)
    diabat_contri = diabat_r + diabat_a + diabat_z
    friction_r = np.nanmean(fric_r, axis=0)
    friction_a = np.nanmean(fric_a, axis=0)
    friction_z = np.nanmean(fric_z, axis=0)
    friction_contri = friction_r + friction_a + friction_z
    horizontal_adv = radial_adv + azimuthal_adv
    total_adv = horizontal_adv + vertical_adv
    total_rhs = total_adv + diabat_contri + friction_contri


    #################################
    #        繪製測試圖
    levels = np.arange(-1,1.01, 0.1)*1e-8
    
    def subplots(ax, var, title, levels):
        c = ax.contourf(cyl_curr.r/1000, cyl_curr.z/1000, var,
                        levels=levels, cmap='bwr', extend='both')
        fig.colorbar(c, ax=ax)
        ax.set_title(title)
    
    fig, ax = plt.subplots(4, 4, figsize=(14,8),dpi=150)
    varlist = [total_lhs, total_rhs, total_lhs - total_rhs, total_rhs*0, 
               radial_adv, azimuthal_adv, vertical_adv, total_adv,
               diabat_r, diabat_a, diabat_z, diabat_contri,  
               friction_r, friction_a, friction_z, friction_contri]
               
    titlelist = ['total LHS', 'total RHS', 'LHS - RHS', '', 
                 'Radial ADV', 'Azimuthal ADV', 'Vertical ADV', 'Total ADV',
                 'Radial DIA', 'Azimuthal DIA', 'Vertical DIA', 'Total DIA',
                 'Radial FRIC', 'Azimuthal FRIC', 'Vertical FRIC', 'Total FRIC'
                 ]
    
    ax = ax.reshape(-1)
    for i in range(16):
        subplots(ax[i], varlist[i], titlelist[i], levels)
    ax = ax.reshape(4,4)
    plt.tight_layout()
    plt.savefig(f'PV_budget{total_t}_{output_time}.png', dpi=120)
    plt.show()


os.system('python default_settings_PV.py')
#讀取設定檔
with open('./settings/settings_PV.json','r') as f:
    settings = json.load(f)

    
calc_cpus = settings['system']['calc_cpus']
setting = settings['wrf']
dt = setting['dt']
Nz = settings['cylindrical_grid']['Nz']
Nr = settings['cylindrical_grid']['Nr']
total_t = 1
buffer = total_t//2
if __name__ == '__main__':
    start_second = datetime.datetime(settings['wrf']['start_yr'], 
                                     settings['wrf']['start_mo'], 
                                     settings['wrf']['start_dy'], 
                                     settings['wrf']['start_hr'], 
                                     settings['wrf']['start_mi'], 
                                     settings['wrf']['start_sc']).timestamp()
    end_second = datetime.datetime(settings['wrf']['end_yr'], 
                                     settings['wrf']['end_mo'], 
                                     settings['wrf']['end_dy'], 
                                     settings['wrf']['end_hr'], 
                                     settings['wrf']['end_mi'], 
                                     settings['wrf']['end_sc']).timestamp()
    dt = settings['wrf']['dt']
    number_of_frame = int((end_second - start_second)/dt) + 1
    t_list = [t for t in range(0+buffer, number_of_frame-buffer)]
    with Pool(processes=12) as pool:
        results = pool.map(plotting_PV_RZ, t_list)
        #results = pool.map(plotting_PV_RZ, t_list)
    
    
    
    

    

    


