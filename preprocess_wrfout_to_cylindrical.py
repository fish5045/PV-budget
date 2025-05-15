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
import wrf
from meteotools.grids import wrfout_grid, cylindrical_grid
from meteotools.fileprocess import get_time_sec
from multiprocessing import Pool

################
## functions
def start_end_second(settings):
    start_second = get_time_sec(settings['start_yr'], settings['start_mo'], 
                                settings['start_dy'], settings['start_hr'], 
                                settings['start_mi'], settings['start_sc'], 0)
    end_second = get_time_sec(settings['end_yr'], settings['end_mo'], 
                              settings['end_dy'], settings['end_hr'], 
                              settings['end_mi'], settings['end_sc'], 0)
    return start_second, end_second

def shift(prevvar, nextvar):
    for j in range(30):
        for i in range(30):
            if prevvar[200+j, 200+i] == nextvar[200, 200]:
                return (j,i)
            elif prevvar[200+j, 200-i] == nextvar[200, 200]:
                return (j,-i)
            elif prevvar[200-j, 200+i] == nextvar[200, 200]:
                return (j,-i)
            elif prevvar[200-j, 200-i] == nextvar[200, 200]:
                return (-j,-i)

def preprocess_cyl(idx):
    wrfdata = wrfout_grid(settings['wrf'], interp_cpus=calc_cpus)
    cyl = cylindrical_grid(settings['cylindrical_grid'], 
                           interp_cpus=calc_cpus)
    wrfdata.set_offset(dt*idx)
        
    #%%
    
    ####################################################
    ###                                              ###
    ###  第 1 區        讀取WRF資料、準備衍伸變數         ###
    ###                                              ###
    ####################################################
    
    ###########################
    ## 1-1   讀取WRF資料與處理
    
    # 讀取WRF資料
    
    wrfdata.read_data('ua', 'va', 'wa', 'MAPFAC_MX',
                      'MAPFAC_MY', 'F', 'pvo', 'XLONG', 'XLAT',
                      'z', 'pressure', 'tk', 'theta', 
                      'RUBLTEN', 'RVBLTEN' ,'H_DIABATIC', 'QVAPOR')
    
    
    # WRF map scale 修正
    wrfdata.correct_map_scale()
    curr_time = wrfdata.wrfout_time
    
    #單位修正
    wrfdata.pressure *= 100   # hPa -> Pa
    wrfdata.pvo *= 1e-6   # PVU -> SI
    
  
    
    #%%
    ####################################################
    ###                                              ###
    ###  第 2 區           內插置圓柱座標               ###
    ###                                              ###
    ####################################################

    ####################################
    ## 2-1   wrfdata內插到直角坐標cyl(z)
    cyl.set_horizontal_location(
        wrfdata.center_xloc, wrfdata.center_yloc)
    
    
    cyl.set_data(wrfdata.dx, wrfdata.dy, wrfdata.z,ver_interp_order=1,
                      u=wrfdata.ua,v=wrfdata.va,w=wrfdata.wa,
                      f=wrfdata.F, p=wrfdata.pressure,
                      T=wrfdata.tk, Th=wrfdata.theta, qv=wrfdata.QVAPOR,
                      fric_u=wrfdata.RUBLTEN, fric_v=wrfdata.RVBLTEN, 
                      LH=wrfdata.H_DIABATIC, PV=wrfdata.pvo)
    

    
    #%%
    ####################################
    ## 2-1   特別處理dPV/dt
    wrfdata_prev = wrfout_grid(settings['wrf'], interp_cpus=calc_cpus)
    wrfdata_prev.set_offset(dt*(idx-1))
    wrfdata_prev.read_data('pvo', 'z', 'XLONG', 'XLAT')
    wrfdata_prev.pvo *= 1e-6   # hPa -> Pa
    
    wrfdata_prev2 = wrfout_grid(settings['wrf'], interp_cpus=calc_cpus)
    wrfdata_prev2.set_offset(dt*(idx-2))
    wrfdata_prev2.read_data('pvo', 'z', 'XLONG', 'XLAT')
    wrfdata_prev2.pvo *= 1e-6   # hPa -> Pa
    
    wrfdata_next = wrfout_grid(settings['wrf'], interp_cpus=calc_cpus)
    wrfdata_next.set_offset(dt*(idx+1))
    wrfdata_next.read_data('pvo', 'z', 'XLONG', 'XLAT')
    wrfdata_next.pvo *= 1e-6   # hPa -> Pa
    
    wrfdata_next2 = wrfout_grid(settings['wrf'], interp_cpus=calc_cpus)
    wrfdata_next2.set_offset(dt*(idx+2))
    wrfdata_next2.read_data('pvo', 'z', 'XLONG', 'XLAT')
    wrfdata_next2.pvo *= 1e-6   # hPa -> Pa
    
                
    yidx, xidx = shift(wrfdata_prev.XLONG, wrfdata.XLONG)
    cyl_prev = cylindrical_grid(settings['cylindrical_grid'], 
                           interp_cpus=calc_cpus)
    cyl_prev.set_horizontal_location(
        wrfdata_prev.center_xloc + xidx*wrfdata_prev.dx, 
        wrfdata_prev.center_yloc + yidx*wrfdata_prev.dy)
    cyl_prev.set_data(wrfdata_prev.dx, wrfdata_prev.dy, wrfdata_prev.z, ver_interp_order=1,
                 PV=wrfdata_prev.pvo)
    
    yidx, xidx = shift(wrfdata_prev2.XLONG, wrfdata.XLONG)
    cyl_prev2 = cylindrical_grid(settings['cylindrical_grid'], 
                           interp_cpus=calc_cpus)
    cyl_prev2.set_horizontal_location(
        wrfdata_prev2.center_xloc + xidx*wrfdata_prev2.dx, 
        wrfdata_prev2.center_yloc + yidx*wrfdata_prev2.dy)
    cyl_prev2.set_data(wrfdata_prev2.dx, wrfdata_prev2.dy, wrfdata_prev2.z, ver_interp_order=1,
                 PV=wrfdata_prev2.pvo)
    
    yidx, xidx = shift(wrfdata.XLONG, wrfdata_next.XLONG)
    cyl_next = cylindrical_grid(settings['cylindrical_grid'], 
                           interp_cpus=calc_cpus)
    cyl_next.set_horizontal_location(
        wrfdata_next.center_xloc - xidx*wrfdata_next.dx, 
        wrfdata_next.center_yloc - yidx*wrfdata_next.dy)
    cyl_next.set_data(wrfdata_next.dx, wrfdata_next.dy, wrfdata_next.z, ver_interp_order=1,
                 PV=wrfdata_next.pvo)
    
    yidx, xidx = shift(wrfdata.XLONG, wrfdata_next2.XLONG)
    cyl_next2 = cylindrical_grid(settings['cylindrical_grid'], 
                           interp_cpus=calc_cpus)
    cyl_next2.set_horizontal_location(
        wrfdata_next2.center_xloc - xidx*wrfdata_next2.dx, 
        wrfdata_next2.center_yloc - yidx*wrfdata_next2.dy)
    cyl_next2.set_data(wrfdata_next2.dx, wrfdata_next2.dy, wrfdata_next2.z, ver_interp_order=1,
                 PV=wrfdata_next2.pvo)
    
    P4 = cyl_next2.PV
    P3 = cyl_next.PV
    P2 = cyl_prev.PV
    P1 = cyl_prev2.PV
    
    cyl.dPVdt = (-P4/12 + P3*2/3 - P2*2/3 + P1/12)/dt
    #%%
    
    ####################################################
    ###                                              ###
    ###  第 3 區               輸出                   ###
    ###                                              ###
    ####################################################

    ####################################
    ## 3-1   wrfdata內插到圓柱坐標cyl(z)
    
    out_file_name = f'./data/PV_pre_cyl_{curr_time}.nc'
    cyl.create_ncfile(out_file_name, wrfdata.time)
    cyl.ncfile_addvar('u', ['vertical', 'tangential', 'radial'], 
                      unit='m/s', description='WE wind')
    cyl.ncfile_addvar('v', ['vertical', 'tangential', 'radial'], 
                      unit='m/s', description='NS wind')
    cyl.ncfile_addvar('w', ['vertical', 'tangential', 'radial'], 
                      unit='m/s', description='vertical velocity')
    cyl.ncfile_addvar('f', ['tangential', 'radial'], 
                      unit='s-1', description='Coriolis parameters')
    cyl.ncfile_addvar('p', ['vertical', 'tangential', 'radial'], 
                      unit='Pa', description='pressure')
    cyl.ncfile_addvar('T', ['vertical', 'tangential', 'radial'], 
                      unit='K', description='temperature')
    cyl.ncfile_addvar('Th', ['vertical', 'tangential', 'radial'], 
                      unit='K', description='potential temperature')
    cyl.ncfile_addvar('qv', ['vertical', 'tangential', 'radial'], 
                      unit='kg/kg', description='mixing ratio of vapor')
    cyl.ncfile_addvar('fric_u', ['vertical', 'tangential', 'radial'], 
                      unit='m/s2', description='WE PBL friction')
    cyl.ncfile_addvar('fric_v', ['vertical', 'tangential', 'radial'], 
                      unit='m/s2', description='NS PBL friction')
    cyl.ncfile_addvar('LH', ['vertical', 'tangential', 'radial'], 
                      unit='K/s', description='latent heating')
    cyl.ncfile_addvar('PV', ['vertical', 'tangential', 'radial'], 
                      unit='m2K/kgs', description='Potential vorticity')
    cyl.ncfile_addvar('dPVdt', ['vertical', 'tangential', 'radial'], 
                      unit='m2K/kgs2', description='Potential vorticity tendency')

    size = os.path.getsize(out_file_name)
    print(f'Complete to preprocess {curr_time}.', flush=True)
    print(f'File: {out_file_name}', flush=True)
    print(f'Size: {size/1024:.1f} KB', flush=True)
    
    from meteotools.calc import FD2
    Th = np.array(wrf.interplevel(wrfdata.theta, wrfdata.z, np.arange(60,6061, 500.)))
    fric_u = np.array(wrf.interplevel(wrfdata.RUBLTEN, wrfdata.z, np.arange(60,6061, 500.)))
    fric_v = np.array(wrf.interplevel(wrfdata.RVBLTEN, wrfdata.z, np.arange(60,6061, 500.)))
    dThdx = FD2(Th, 1000., 2)
    dThdy = FD2(Th, 1000., 1)
    dThdz = FD2(Th, 500., 0)
    cross_x = -FD2(fric_v, 500, 0)
    cross_y = FD2(fric_u, 500, 0)
    cross_z = FD2(fric_v, 1000, 2) - FD2(fric_u, 1000, 1)
    friction_x = cross_x*dThdx
    friction_y = cross_y*dThdy
    friction_z = cross_z*dThdz
    
    
    


os.system('python default_settings_PV_pre.py')
#讀取設定檔
with open('./settings/settings_PV_pre.json','r') as f:
    settings = json.load(f)

    
calc_cpus = settings['system']['calc_cpus']
dt = settings['wrf']['dt']
wrf.omp_set_num_threads(calc_cpus)
start_second, end_second = start_end_second(settings['wrf'])
max_files = int((end_second - start_second)/dt) + 1

if __name__ == '__main__':
    idx_list = [i for i in range(max_files)]
    with Pool(processes=12) as pool:
        results = pool.map(preprocess_cyl, idx_list)
    #preprocess_cyl(idx_list[0])


