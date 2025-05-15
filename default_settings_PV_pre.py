import json
import numpy as np

d = dict()
d['wrf'] = dict()
d['cylindrical_grid'] = dict()
d['cartesian_grid'] = dict()
d['system'] = dict()


dwrf = d['wrf']                #wrfout 設定
dwrf['wrfout_dir'] = r'C:\Users\fish5045\Desktop\work\wrfdata'      #wrfout 位置
dwrf['start_yr'] = 2019        #wrfout 診斷起始年
dwrf['start_mo'] = 8           #              月
dwrf['start_dy'] = 8           #              日
dwrf['start_hr'] = 11          #              十
dwrf['start_mi'] = 00           #              分
dwrf['start_sc'] = 0           #              秒
dwrf['end_yr'] = 2019          #           終止年
dwrf['end_mo'] = 8             #              月
dwrf['end_dy'] = 8             #              日
dwrf['end_hr'] = 14            #              時
dwrf['end_mi'] = 00             #              分
dwrf['end_sc'] = 0             #              秒
dwrf['domain'] = 4             #              domain編號
dwrf['dt']     = 300                          #2wrfout的時間間格dt

#直角坐標系網格設定 (用於疊代擾動氣壓)
dcyl = d['cylindrical_grid']              #圓柱座標設定
dcyl['Nr'] = 151                          #徑向格數
dcyl['Ntheta'] = 360                       #切向格數
dcyl['Nz'] = 40                           #垂直格數
dcyl['dr'] = 1000                         #徑向網格大小 (m)
dcyl['dtheta'] = 2*np.pi/dcyl['Ntheta']   #切向網格大小 (rad)
dcyl['dz'] = 500                          #垂直網格大小 (m)
dcyl['z_start'] = 60                      #垂直起始位置 (m)
dcyl['r_start'] = 0                       #徑向起始位置 (m)
dcyl['theta_start'] = 0                   #切向起始位置 (rad)
dcyl['vertical_coord_type'] = 'z'         #垂直座標
dcyl['dt'] = 300                          #tendency 使用的dt



#系統相關設定
dsys = d['system']
dsys['calc_cpus'] = 2
dsys['interp_cpus'] = 8


with open('./settings/settings_PV_pre.json', 'w') as f:
    #json.dump(d, f)
    f.write(json.dumps(d,indent=4))

'''
with open('settingss.json', 'r') as f:
    k = json.load(f)
'''


