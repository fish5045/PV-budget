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
import time as tm
import matplotlib as mpl
import os
import datetime
from meteotools.grids import cylindrical_grid
from meteotools.fileprocess import get_wrfout_time, mkdir
from meteotools.calc import FD2
from multiprocessing import Pool

################
## functions
def prev_curr_next_time(settings,shift=0):
    prev_time2 = get_wrfout_time(settings['start_yr'], settings['start_mo'], 
                                settings['start_dy'], settings['start_hr'], 
                                settings['start_mi'], settings['start_sc'],
                                -settings['dt']*2+shift)
    prev_time = get_wrfout_time(settings['start_yr'], settings['start_mo'], 
                                settings['start_dy'], settings['start_hr'], 
                                settings['start_mi'], settings['start_sc'],
                                -settings['dt']+shift)
    curr_time = get_wrfout_time(settings['start_yr'], settings['start_mo'], 
                                settings['start_dy'], settings['start_hr'], 
                                settings['start_mi'], settings['start_sc'],
                                0+shift)
    next_time = get_wrfout_time(settings['start_yr'], settings['start_mo'], 
                                settings['start_dy'], settings['start_hr'], 
                                settings['start_mi'], settings['start_sc'],
                                settings['dt']+shift)
    next_time2 = get_wrfout_time(settings['start_yr'], settings['start_mo'], 
                                settings['start_dy'], settings['start_hr'], 
                                settings['start_mi'], settings['start_sc'],
                                settings['dt']*2+shift)
    return prev_time2, prev_time, curr_time, next_time, next_time2

def compute_PV_budget(curr_time):
    ####################################################
    ###                                              ###
    ###  第 1 區    讀取圓柱座標資料、準備衍伸變數        ###
    ###                                              ###
    ####################################################
    
    ###########################
    ## 1-1   讀取圓柱座標資料
    
    #讀取設定檔、初始化網格資料
    print(f'{curr_time}')
    
    cyl_curr = cylindrical_grid(f'./data/PV_pre_cyl_{curr_time}.nc')
    cyl_curr.ncfile_readvars('u', 'v', 'w', 'f', 'p', 'T', 'Th', 'LH', 'dPVdt',
                             'PV', 'fric_u', 'fric_v', 'qv')
    
    ###########################
    ## 1-2   衍伸變數
    cyl_curr.set_horizontal_location(0,0)
    
    cyl_curr.calc_fric_vr_vt()

    cyl_curr.calc_abs_vorticity_3D()
    cyl_curr.calc_rho()
    cyl_curr.calc_density_potential_temperature()
    
    #%%    
    #%%
    ####################################################
    ###                                              ###
    ###  第 2 區        PV收支分析                     ###
    ###                                              ###
    ####################################################

    ####################################
    ## 2-1   平均場、擾動場 (假設沒有Fz、不考慮水物終端速度)
    
    
    dr = cyl_curr.dr
    dtheta = cyl_curr.dtheta
    dz = cyl_curr.dz
    r = cyl_curr.r.reshape(1,1,-1)
    z = cyl_curr.z
    PV = cyl_curr.PV
    vt = cyl_curr.vt
    vr = cyl_curr.vr
    w = cyl_curr.w
    LH = cyl_curr.LH
    rho = cyl_curr.rho
    zeta_r = cyl_curr.zeta_r
    zeta_theta = cyl_curr.zeta_t
    zeta_z = cyl_curr.abs_zeta_z
    Th = cyl_curr.Th
    fric_vr = cyl_curr.fric_vr
    fric_vt = cyl_curr.fric_vt
    #print(np.nanmax(fric_vt), np.nanmax(fric_vr))

    
    ####################################
    ## 2-2   RHS (假設沒有Fz、不考慮水物終端速度)
    # 注意微分維度 var (z, theta, r)
    #            sym_var (z, r)
    #            asy_var (z, theta_r)
    
    
    def gradient(var, dr, dtheta, dz, r):
        gradient_r = FD2(var, dr, 2)
        gradient_theta = FD2(var, dtheta, 1)/r.reshape(1,1,-1)
        gradient_z = FD2(var, dz, 0)
        return gradient_z, gradient_theta, gradient_r
    
    def cross(ar, aa, az, br, ba, bz):
        cross_r = aa*bz - az*ba
        cross_a = az*br - ar*bz
        cross_z = ar*ba - aa*br
        return cross_z, cross_a, cross_r
        
    
    dThdz, dThdtheta, dThdr = gradient(Th, dr, dtheta, dz, r)
    dPVdz, dPVdtheta, dPVdr = gradient(PV, dr, dtheta, dz, r)
    dLHdz, dLHdtheta, dLHdr = gradient(LH, dr, dtheta, dz, r)
    drhodz, drhodtheta, drhodr = gradient(rho, dr, dtheta, dz, r)
    
    
    cross_r = -FD2(fric_vt, dz, 0)
    cross_a = FD2(fric_vr, dz, 0)
    cross_z = (FD2(r*fric_vt, dr, 2) - FD2(fric_vr, dtheta, 1))/r
    
    adv_a = -vt*dPVdtheta
    adv_r = -vr*dPVdr
    adv_z = -w*dPVdz
    dia_r = zeta_r*dLHdr/rho
    dia_a = zeta_theta*dLHdtheta/rho
    dia_z = zeta_z*dLHdz/rho
    fric_r = dThdr*cross_r/rho
    fric_a = dThdtheta*cross_a/rho
    fric_z = dThdz*cross_z/rho
    
    ####################################
    ## 2-3   alias 變數別名
    
    cyl_curr.lhs = cyl_curr.dPVdt
    cyl_curr.adv_r = adv_r
    cyl_curr.adv_a = adv_a
    cyl_curr.adv_z = adv_z
    cyl_curr.dia_r = dia_r
    cyl_curr.dia_a = dia_a
    cyl_curr.dia_z = dia_z
    cyl_curr.fric_r = fric_r
    cyl_curr.fric_a = fric_a
    cyl_curr.fric_z = fric_z

    mkdir(f'./output/')
    out_file_name = f'./output/PV_budget_cyl_{curr_time}.nc'   #滑動平均最後一個時間點
    cyl_curr.create_ncfile(out_file_name)
    cyl_curr.ncfile_addvar('lhs', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Local PV Tendency')
    cyl_curr.ncfile_addvar('adv_r', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgss', description='Radial PV Advection')
    cyl_curr.ncfile_addvar('adv_a', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Azimuthal PV Advection')
    cyl_curr.ncfile_addvar('adv_z', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Vertical PV Advection')
    cyl_curr.ncfile_addvar('dia_r', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Radial Diabat Contriction')
    cyl_curr.ncfile_addvar('dia_a', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Azimuthal Diabat Contriction')
    cyl_curr.ncfile_addvar('dia_z', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Vertical Diabat Contriction')
    cyl_curr.ncfile_addvar('fric_r', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Radial Friction Contriction')
    cyl_curr.ncfile_addvar('fric_a', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Azimuthal Friction Contriction')
    cyl_curr.ncfile_addvar('fric_z', ['vertical', 'tangential', 'radial'], 
                      unit='Km2/kgs', description='Vertical Friction Contriction')
    
    
    return 



os.system('python default_settings_PV.py')
#讀取設定檔
with open('./settings/settings_PV.json','r') as f:
    settings = json.load(f)

    
calc_cpus = settings['system']['calc_cpus']
dt = settings['wrf']['dt']
prev_time2, prev_time, curr_time, next_time, next_time2 = \
    prev_curr_next_time(settings['wrf'])

if __name__ == '__main__':
    total_t = 1
    cyl_curr = cylindrical_grid(f'./data/PV_pre_cyl_{curr_time}.nc')
    cyl_curr.ncfile_readvars('u', 'v', 'w')
    cyl_curr.calc_vr_vt()
    Nz = cyl_curr.Nz
    Nr = cyl_curr.Nr
    Ntheta = cyl_curr.Ntheta
    
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
    curr_time_list = []
    for t in range(0, number_of_frame):
        prev_time2, prev_time, curr_time, next_time, next_time2 = \
            prev_curr_next_time(settings['wrf'],shift=t*dt)
        curr_time_list.append(curr_time)
    
    with Pool(processes=12) as pool:
        results = pool.map(compute_PV_budget, curr_time_list)
    #compute_PV_budget(curr_time_list[0])

    

    


