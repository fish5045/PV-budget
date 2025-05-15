# 無標題

# 安裝

### 環境檢查

使用前請確保您的環境已安裝以下套件:

numpy scipy wrf-python netCDF4 json (opencv)

| numpy | scipy | wrf-python | netCDF4 | python |
| --- | --- | --- | --- | --- |

註: opencv可選用安裝，若未安裝則無法使用meteotools.pictureprocess.movie_make。

使用前請確保您的系統已安裝以下軟體:

gfortran [https://www.mingw-w64.org/downloads/#mingw-builds](https://www.mingw-w64.org/downloads/#mingw-builds)

將gfortran加入PATH，使在指令視窗下可直接輸入gfortran即找到此編譯器。

請參考: [安裝編譯器並將編譯器加入系統路徑(PATH)](https://www.notion.so/PATH-1f16f2851fb18054b4a5c8f0ada1aee1?pvs=21)

本套件是在以下環境測試:

| 名稱 | 版本 |
| --- | --- |
| python | 3.8.12 |
| numpy | 1.24.4 |
| scipy | 1.9.1 |
| nctCDF4 | 1.6.4 |
| wrf-python | 1.3.4.1 |
| gfortran 64-bit | 14.2.0 |
| Windows 11 64-bit |  |

# 建置

下載PV_budgets後，請先至meteotools資料夾內執行一次python build.py，將依照您的作業系統編譯套件需要的函式庫。

請確保gfortran已安裝至您的系統，並將gfortran加入PATH，並可直接由指令呼叫。

# 使用

1. 編輯default_settings_PV_pre.py，進行WRFOUT的位置、起訖時間點、資料時間間格設定，及目標的圓柱座標網格系統設定。
2. 編輯default_settings_PV.py，將裡面的設定設成default_settings_PV_pre.py的內容。這裡的起迄時間點為診斷的起迄時間，時段必須小於等於default_settings_PV_pre.py所設定的時段。
3. 執行preprocess_wrfout_to_cylindrical.py。
4. 執行main_qPV.py。
5. 執行PV_plotting.py。

### 註:

1. default_settings_PV_pre.py、default_settings_PV.py內可設定calc_cpus與，用於計算同時處理的資料筆數、由WRFOUT轉換至圓柱座標內插時所用的核心數，在硬體許可下以縮短等待時間。
建議calc_cpus最大值不要超過實體核心數。
請考量您電腦的記憶體大小，測試時將WRFOUT(36085500格)轉換為圓柱座標(4348800格)所需的記憶體為2.5GB*calc_cpus。
如果您的WRFOUT資料存在傳統硬碟，不建議將calc_cpus設定超過4。
2. PV_plotting.py內可以設定繪圖時取滑動平均的時段，請修改變數total_t，當total_t設為7時，將使用前3筆、後3筆、與本次時間的7筆資料進行平均。