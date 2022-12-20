import numpy as np
import pandas as pd
from scipy import interpolate
#import mod_fvcom
#import sys
def interp(df,siglay,date_ary,var_name):
    #out_df = pd.DataFrame()      
    df = df.dropna(subset=['tp','sl'])      #nanを落とす
    #levmax = max(df['lev'])
    out_df = pd.DataFrame()
    depth = [df[df['datetime'] == date].tail(1).iat[0,7]  for date in date_ary.values if len(df[df['datetime'] == date])!=0]
    mean_depth = sum(depth)/len(depth)
    #with open('./data/'+place+var_name+'depth.txt','w') as f:
    #    print(depth,f)
    for date in date_ary.values:
        
        #mask = (df['datetime'] == date)
        tmp = df[df['datetime'] == date]
        if var_name == 'salinity':
            tmp =tmp[tmp['sl'] >=15]
        else:
            tmp = tmp[tmp['tp']<=40]
        if len(tmp) == 0:
            l= [np.nan*len(siglay)]
        else:
            x  = np.array(tmp['depth'])
            if tmp.tail(1).iat[0,7]>mean_depth/2 and len(tmp)>1: #最大水深が浅すぎる場合は除く
                    if var_name == 'temperature':
                        y  = np.array(tmp['tp'])
                    else:
                        y  = np.array(tmp['sl'])
                    #print(x,y)
                    #sigma_d = np.array([max(x)*siglay[i] for i in range(len(siglay))])
                    sigma_d = np.array([mean_depth*siglay[i] for i in range(len(siglay))])
                    #sigma_d = np.array(sigma_d)

                    f = interpolate.interp1d(x, y,kind='linear',axis=-1,copy=True,bounds_error=None, \
                                            fill_value='extrapolate',assume_sorted=False)
                    var_new = f(sigma_d)
                    l= var_new
            else:
                l= [np.nan*len(siglay)]
        l = pd.Series(l,name=date)

        out_df = pd.concat([out_df,l],axis=1)
       # out_df = out_df.T
    return out_df

date_index = pd.date_range("2020-01-01 00:00:00", periods=24*366, freq="H") #specify start date and output frequency
date_ary = date_index.to_series().dt.strftime("%Y-%m-%d %H:%M:%S")
obcfile = '../input_testcase_interp/TokyoBay_obc.dat'
#f = '../run442/Tokyo_0001.nc'
#fvcom = mod_fvcom.Fvcom2D(ncfile_path=f, obcfile_path=obcfile, m_to_km=True, offset=False)

print('lets begin!')


#f = '/work/gy29/y29007/Github/data/csv/Mpos_S_' + key +'_2020.csv' #OBCX
variables = ['salinity','temperature']
variables = ['temperature']
stations = ['chiba1buoy','chibaharo','kawasaki','urayasu']
stations = ['chiba1buoy','chibaharo']
for variable in variables:
    for place in stations:
        f = './data/Mpos_S_'+place+'_2020.csv'
        df = pd.read_csv(f)

        siglay = [0.00222222222222222,0.00888888888888889,0.0200000000000000,0.0355555555555556,0.0555555555555556,0.0800000000000000,0.108888888888889,0.142222222222222,0.180000000000000,0.222222222222222,0.268888888888889,0.320000000000000,0.375555555555556,0.435555555555556,0.500000000000000,0.564444444444444,0.624444444444444,0.680000000000000,0.731111111111111,0.777777777777778,0.820000000000000,0.857777777777778,0.891111111111111,0.920000000000000,0.944444444444444,0.964444444444444,0.980000000000000,0.991111111111111,0.997777777777778]
        #siglay = -fvcom.siglay[:,index].data
        #print(siglay)
        regress = interp(df,siglay,date_ary,variable)
        regress.to_csv('./data/interp/'+place+'_'+variable+'sigma2_2020.csv')
        print(f"done.{place}, {variable}")
"""
siglay = [0.00222222222222222,0.00888888888888889,0.0200000000000000,0.0355555555555556,0.0555555555555556,0.0800000000000000,0.108888888888889,0.142222222222222,0.180000000000000,0.222222222222222,0.268888888888889,0.320000000000000,0.375555555555556,0.435555555555556,0.500000000000000,0.564444444444444,0.624444444444444,0.680000000000000,0.731111111111111,0.777777777777778,0.820000000000000,0.857777777777778,0.891111111111111,0.920000000000000,0.944444444444444,0.964444444444444,0.980000000000000,0.991111111111111,0.997777777777778]

        siglay = [0.01666667 ,0.05      , 0.08333334 ,0.11666667 ,0.15000001 ,0.18333334, \
        0.21666668 ,0.25000001, 0.28333335 ,0.31666668 ,0.35000001 ,0.38333336, \
        0.4166667  ,0.45000003, 0.48333335 ,0.51666668 ,0.55000004 ,0.58333337, \
        0.6166667  ,0.65000004, 0.68333337 ,0.7166667  ,0.75000003 ,0.78333339, \
        0.81666672 ,0.85000005, 0.88333338 ,0.91666672 ,0.95000005 ,0.98333335]
"""