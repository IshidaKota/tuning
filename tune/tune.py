
"""
Usage:
python tune.py run_directory input.file variable_name
line 474~ user difinition 
without extention
vaiable_name:salinity or temperature
"""

import pandas as pd
import numpy as np
import math
import sys
import os
#from functools import partial, partialmethod
from pyproj import Proj
### import xarray as xr
#import netCDF4
import matplotlib.pyplot as plt
import holoviews as hv
#import subprocess
from scipy import interpolate
#hv.extension('bokeh', 'matplotlib')
import mod_fvcom
from matplotlib.colors import Normalize
from PyFVCOM.plot import Depth#,Time
#from PyFVCOM.grid import unstructured_grid_depths
from PyFVCOM.read import FileReader
from holoviews import opts
import matplotlib.animation as animation
import datetime
from jdcal import gcal2jd, jd2gcal ,MJD_0
from matplotlib import dates as mdates
def find_closest_node(stn):
    """
    moniteringpostから最近傍の場所を抽出
    """
    index = {}
    p = Proj(proj='utm',zone=54,ellps='WGS84', preserve_units=False)
    
    for key in stn.keys():
        st_m = p(stn[key][0],stn[key][1]) #lon,lat→utm(m)
        st_km= st_m[0]/1000,st_m[1]/1000 #utm(m)→utm(km)
        dist = [math.sqrt((fvcom.x[i]-st_km[0])**2 +(fvcom.y[i]-st_km[1])**2) for i in range(len(fvcom.x))] #全nodeからの距離を計算、listで保存
        #listの中から最小値を見つける
        res ,min_dist= dist.index(min(dist)),min(dist) #indexを取得
        index[key] = res  
        print(f"Node {res+1} is {min_dist}km to the station:{key},station is at {stn[key]}")
         ##+1 :SMS上でnodeは1始まり、listのindexは0始まり
        
    return index

def calc_corr(x,y):
    s1 = pd.Series(x)
    s2 = pd.Series(y)
   # print(s1,s2)
    res = s2.corr(s1)
    return res
    



def interp(df,siglay,date_ary,var_name):
    #out_df = pd.DataFrame()      
    df = df.dropna(subset=['tp','sl'])      #nanを落とす
    levmax = max(df['lev'])
    interp0 = {}
    for date in date_ary.values:
        #mask = (df['datetime'] == date)
        tmp = df[df['datetime'] == date]
        if len(tmp) == 0:
            interp0[date] = np.zeros(len(siglay))
        else:
            x  = np.array(tmp['depth'])
            if x.size >= levmax/2:
                if var_name == 'temperature':
                    y  = np.array(tmp['tp'])
                else:
                    y  = np.array(tmp['sl'])
                    
                sigma_d = np.array([max(x)*siglay[i] for i in range(len(siglay))])
                #sigma_d = np.array(sigma_d)

                f = interpolate.interp1d(x, y,kind='linear',axis=-1,copy=True,bounds_error=None, \
                                           fill_value='extrapolate',assume_sorted=False)
                var_new = f(sigma_d)
                interp0[date] = var_new
            else:
                interp0[date] = np.zeros(len(siglay))
    return interp0
    
   


def plotter(stn_node,date_ary,var_name): 
    
    if var_name == 'temperature':
        var = fvcom.variables['temp']
    elif var_name == 'salinity':
        var = fvcom.variables['salinity']
    else:
        raise ValueError(f"invalid variable '{var_name}'")

    corr = {};rmse = {};df = pd.DataFrame()
    for key in stn_node:
        index = stn_node[key]
        #read Mpos_S
        #f = '/work/gy29/y29007/Github/data/csv/Mpos_S_' + key +'_2020.csv'
        #df = pd.read_csv(f,usecols=['depth','tp','sl','datetime','lev'])
        regress = regress_dict_tp[key] #dictのdict形式でまとまっている
        sigma =[2,15,len(siglay)-2]
        for k in sigma:
            res = np.zeros(len(date_ary))
            i = 0
            for date in date_ary.values:
                if regress[date][0] == 0:               
                    res[i] = np.nan
                else: 
                    res[i] =regress[date][k]
                #debug
                if res[i]>50: 
                    print(f"{date} sal or temp value is {res[i]}")
                    res[i] = np.nan
                    
                if  res[i]<10:
                 #   print(f"{date} sal or temp value is {res[i]}")
                    res[i] = np.nan
                    

                i+=1
            head = key+str(k)
            df[head] = res
            
    #plot
            fig = plt.figure(figsize=(6,3))
            ax = fig.add_subplot(1,1,1,xlabel = 'month(2020)',ylabel='{}'.format(var_name))
            ax.plot(var[:,k, index],linewidth=2,label='simulation') #float64 temp(time, siglay, node)
            ax.plot(res[:],linewidth=2,label='observation')
            fig.suptitle('{}@{},sigma = {}/30'.format(var_name,key,k))
            ax.legend()
            plt.grid()
            if var_name == 'temperature':
                plt.ylim(-5,30)
            else :
                plt.ylim(10,35)
            #plt.xticks([24*30*i for i in range(12)], [30*i for i in range(12)])
            plt.xticks([30*i for i in range(12)], [i for i in range(1,13)])
            plt.plot()
            os.makedirs('./png/'+filename,exist_ok=True)
            fig.savefig('./png/'+filename+'/'+var_name+'_'+key+'_'+str(k)+'_'+filename+'.png',\
               bbox_inches='tight', pad_inches=0.1,dpi=600) 
            print('save figure@sigma={},stn={}'.format(k,key))
            
    #calc RMSE and CC
            #RMSE
            print('start calculate RMSE and CORR')
            X = []
            Y = []
            diff = 0
            for i in range(min(len(fvcom.time),len(res))):
                #print(np.isnan(res[i]))
                if (np.isnan(res[i]) is False) or (res[i]>0):
                    #print(res[i])
                    X.append(var[i,k,index].item())
                    Y.append(res[i])
                    diff += (var[i,k,index].item()-res[i])**2
            #print(f"Var is {diff}") #debug
            Rmse = math.sqrt(diff/len(Y)) #root 1/n

            #Corr
            #print(f"lenx={len(X)},leny={len(Y)},{X[:10]},{Y[:10]}")
            Corr = calc_corr(X,Y)

            corr[key] = Corr
            rmse[key] = Rmse
    #os.makedirs('./result',exist_ok=True)
    df.to_csv('./df/w15r15'+filename+var_name+'_simulation.csv')
    return rmse,corr

def plotter2(var_name,place,index):
    #which coordinate does sigma use?

    if sigmoid:
        f = './data/interp/'+place+'_'+var_name+'sigma2_2020.csv'
    else:
        f = './data/interp/'+place+'_'+var_name+'_2020.csv'
    df = pd.read_csv(f)
    #make time series
    dates_float = [jd2gcal(MJD_0,fvcom.Itime[i]) for i in range(len(fvcom.Itime))] #ymd
    dates = [datetime.datetime(dates_float[i][0],dates_float[i][1],dates_float[i][2],\
        hour=int(fvcom.Itime2[i]/3600000))for i in range(len(fvcom.Itime))] #ymd+h
    column = pd.Series(list(df.columns[1:]))
    column = pd.to_datetime(column)
    
    if var_name == 'temperature':
        var = fvcom.variables['temp']
    elif var_name == 'salinity':
        var = fvcom.variables['salinity']
    corr = [];rmse = [];df_sim = pd.DataFrame()

    sigma =[2,15,28]
    layer = ['upper','middle','bottom']
    fig,axs = plt.subplots(1,3,sharex='all',figsize=(20,4))
    kokyo_df = pd.read_csv('./data/kokyo/'+place+'.csv')
    kokyo_df['date2'] =[ datetime.datetime.strptime(kokyo_df['date'][i],'%Y-%m-%d %H:%M:%S') for i in range(len(kokyo_df))]
    #print(kokyo_df.head())
    for i, k in enumerate(sigma):
        res = np.array(df.iloc[k,:])
       # print(type(column[0]),type(dates[0]),type(kokyo_df['date'][0]))
        #plot
        if variable == 'salinity':
            if layer[i] == 'upper':
                axs[i].plot(kokyo_df['date2'][:18:2],kokyo_df['salinity'][:18:2],'ro',label='kokyo')
            elif layer[i] == 'bottom':
                axs[i].plot(kokyo_df['date2'][1:18:2],kokyo_df['salinity'][1:18:2],'ro',label='kokyo')
  
        axs[i].plot(dates,var[:,k, index],linewidth=2,label='simulation') 
        axs[i].plot(column,res[1:],linewidth=2,label='observation')

        axs[i].plot
        axs[i].xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))

        axs[i].set_title('{} of {}, {}'.format(var_name,place,layer[i]))
        axs[i].legend()
        axs[i].grid()
        if var_name == 'salinity':
            axs[i].set_ylim(15,35)
        else:
            axs[i].set_ylim(0,35)

        os.makedirs('./png/'+filename,exist_ok=True)
      
                    
        #calc RMSE and CC
                    #RMSE
        #print('start calculate RMSE and CORR')
        X = []
        Y = []
        diff = 0
        #for i in range(min(len(fvcom.time),len(res))):
        #    #print(np.isnan(res[i]))
        #    if res[i]>0: # Time_list[i] in regress.keys() and 
        #        #print(res[i])
        #        X.append(var[i,k,index].item())
        #        Y.append(res[i])
        #        diff += (var[i,k,index].item()-res[i])**2
        #print(f"Var is {diff}") #debug
        for i, colum in enumerate(column):
            if colum in dates and not(np.isnan(res[i])) :
                j = dates.index(colum)
                X.append(var[j,k,index].item())
                Y.append(res[i])
                diff += (var[j,k,index].item()-res[i])**2
        Rmse = math.sqrt(diff/len(Y)) #root 1/n

        #Corr
        #print(f"lenx={len(X)},leny={len(Y)},{X[:10]},{Y[:10]}")
        Corr = calc_corr(X,Y)
        corr.append([place,k,Corr])
        rmse.append([place,k,Rmse])

        #save
        df_sim['date'] = dates
        df_sim[str(k)] = var[:,k,index]
    df_sim.to_csv('./df/result/'+place+var_name+'_'+filename+'_.csv')
    print('./df/result/'+place+var_name+'_'+filename+'_.csv')
    #print(df_sim.head())
    os.makedirs('./png/'+filename,exist_ok=True)
    fig.savefig('./png/'+filename+'/'+var_name+'_'+place+'.png',\
            bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.close()
 

    return corr,rmse


#################make velocity map#########################

def make_animation(filename):
    print('start making animation of velocity...')
    from datetime import timedelta
    from datetime import datetime as dt

    sigma = [2,15,28]
    for sigma in sigma:
        fig = plt.figure(facecolor="w")
        ax = fig.add_subplot(1, 1, 1, aspect="equal",ylabel ="Y(km)", xlabel="X(km)")
        x = fvcom.verts['x']
        y = fvcom.verts['y']
        tris = fvcom.tris
        cx = [sum(x[tris.iloc[i,:]])/3 for i in range(len(tris))]
        cy = [sum(y[tris.iloc[i,:]])/3 for i in range(len(tris))]

        #sigma = 0
        u = fvcom.variables['u'][:,sigma,:]
        v = fvcom.variables['v'][:,sigma,:]

        q = plt.quiver(cx,cy,u[0,:],v[0,:], color='red',width=0.005,\
                    angles='xy', scale_units='xy', scale=0.1,label='velocity')

        def update_quiver(t, q, cx, cy):
            ut = u[t,:]
            vt = v[t,:]
            q.set_UVC(ut,vt)
            date = dt.strptime("2020/1/1", '%Y/%m/%d')+ timedelta(days=t) 
            plt.title('sigma={},date={}'.format(sigma,date))
            return q,

        plt.plot()
        ani = animation.FuncAnimation(fig, update_quiver,fargs=(q, cx,cy),interval=200, blit=False)
        ani.save("./png/v/velocity_"+str(sigma)+"_"+filename+".gif", dpi =200,writer="imagemagick")
        print(f"sigma={sigma},done.")
        plt.close()
###############################################################
def plot_outer_sea(filename):
    print("let's confirm outer sea")
    index = 774
    if variable == 'salinity':
        var = fvcom.variables['salinity'][:,:,:].data
    else:
        var = fvcom.variables['temp'][:,:,:].data

    fig = plt.figure(figsize=(6,3));ax = fig.add_subplot(1,1,1,xlabel = 'month(2020)',ylabel='{}'.format(variable))
    for k in [2,15,28]:
        ax.plot(var[:,k, index],linewidth=2,label='sigma'+str(k)) #float64 temp(time, siglay, node)
    ax.legend();plt.grid()
    plt.xticks([30*i for i in range(12)], [i for i in range(1,13)])
    os.makedirs('./png/'+filename,exist_ok=True)
    fig.savefig('./png/'+filename+'/'+variable+'_outer.png',\
    bbox_inches='tight', pad_inches=0.1,dpi=300)  
    plt.close()
    print('done.')
###################################
def plot_map(filename):
    stn_node['miura'] =801
    path = './png/map/'
    switch = 1
    if switch == 1:
        print('start plotting map...')
        # station map
        node1=stn_node['chiba1buoy']; node2=stn_node['chibaharo'];node3=stn_node['kawasaki'] ;node4=stn_node['urayasu'];node5=stn_node['miura']
        coords1 = [[fvcom.variables['x'][node1].item(), fvcom.variables['y'][node1].item()]]
        coords2 = [[fvcom.variables['x'][node2].item(), fvcom.variables['y'][node2].item()]]
        coords3 = [[fvcom.variables['x'][node3].item(), fvcom.variables['y'][node3].item()]]
        coords4 = [[fvcom.variables['x'][node4].item(), fvcom.variables['y'][node4].item()]]
        coords5 = [[fvcom.variables['x'][node5].item(), fvcom.variables['y'][node5].item()]]
        #outer 45 484 589 987
        node6 =45;node7=484;node8=589,;node9=987
        l1 =(fvcom.variables['x'][node6].item(), fvcom.variables['y'][node6].item())
        l2 =(fvcom.variables['x'][node7].item(), fvcom.variables['y'][node7].item())
        l3 =(fvcom.variables['x'][node8].item(), fvcom.variables['y'][node8].item())
        l4 =(fvcom.variables['x'][node9].item(), fvcom.variables['y'][node9].item())
        #opts.defaults(opts.Points(tools=['hover']))
        path1 = hv.Path([l1,l2])
        path2 = hv.Path([l2,l3])
        path3 = hv.Path([l3,l4])
        mesh = fvcom.plot_mesh(line_color='blue', line_width=0.01)
        points1 = fvcom.plot_point(coords1, color='#63000f', marker='x', size=12,label='chiba1buoy')
        points2 = fvcom.plot_point(coords2, color='green', marker='x', size=12,label='chibaharo')
        points3 = fvcom.plot_point(coords3, color='red', marker='x', size=12,label='kawasaki')
        points4 = fvcom.plot_point(coords4, color='blue', marker='x', size=12,label='urayasu')
        points5 = fvcom.plot_point(coords5, color='blue', marker='x', size=12,label='miura')
        p = mesh * points1 * points2 * points3 * points4 * points5 * fvcom.plot_coastline(color='black',linewidth=0.01) *path1 *path2 * path3
        p_opts = opts.Points(frame_width=500, aspect='equal', xlim=(340, 460))
        p.opts(p_opts)

        
        hv.save(p, path+'st_location_outer.png', fmt='png',dpi=300, backend='matplotlib')
        print('done.')

#################USE PYFVCOM############################################
#vertical plot
def plot_vertical(filename):
    print('start making vertical contour...')
    f = '../'+run+'/'+filename+ '_0001.nc'
    fvcom_pyfvcom = FileReader(f, variables=['zeta', 'salinity','temp','kh','viscofh'], dims={'time': slice(0,len(fvcom.time))})
    ##### Start set parameters
    timestep_base = [30,90,150,210,240,270];scale_factor = fvcom.time[1]-fvcom.time[0]
    timestep = [int(timestep_base[i]/scale_factor) for i in range(len(timestep_base)) if int(timestep_base[i]*scale_factor)<len(fvcom.time)]
    #print(timestep)
    #outer 45 484 589 987
    #inner 286 2558 2665 ---
    lon_lat1, lon_lat2 = (140.03767, 35.6040883), (139.698287, 35.337028) 
    positions = np.array((lon_lat1, lon_lat2)) 
    indices, distances = fvcom_pyfvcom.horizontal_transect_nodes(positions)
    print(f"start={indices[0]}, end={indices[-1]}")
    #補正
    p = Proj(proj='utm',zone=54,ellps='WGS84', preserve_units=False)
    x1,y1 =p(lon_lat1[0],lon_lat1[1])[0],p(lon_lat1[0],lon_lat1[1])[1]
    x2,y2 =p(lon_lat2[0],lon_lat2[1])[0],p(lon_lat2[0],lon_lat2[1])[1]
    true_dist =math.sqrt((x1-x2)**2 + ((y1-y2)**2))
    for i in range(len(distances)):
        distances[i] = true_dist*distances[i]/distances[-1]

    lon_lat1, lon_lat2 =(139.698287, 35.337028)  , (139.795921, 35.2481378) ## Estuary example
    x1,y1 =p(lon_lat1[0],lon_lat1[1])[0],p(lon_lat1[0],lon_lat1[1])[1]
    x2,y2 =p(lon_lat2[0],lon_lat2[1])[0],p(lon_lat2[0],lon_lat2[1])[1]
    true_dist1 =math.sqrt((x1-x2)**2 + ((y1-y2)**2))
    positions = np.array((lon_lat1, lon_lat2)) 
    indices1, distances1 = fvcom_pyfvcom.horizontal_transect_nodes(positions)
    print(f"start1={indices1[0]}, end1={indices1[-1]}")
    for i in range(len(distances1)):
        distances1[i] = (true_dist1*distances1[i]/distances1[-1]) +true_dist
    



    if outer:


        lon_lat1, lon_lat2 =(139.795921, 35.2481378) , (139.574171, 34.853014) 
        x1,y1 =p(lon_lat1[0],lon_lat1[1])[0],p(lon_lat1[0],lon_lat1[1])[1]
        x2,y2 =p(lon_lat2[0],lon_lat2[1])[0],p(lon_lat2[0],lon_lat2[1])[1]
        true_dist2 =math.sqrt((x1-x2)**2 + ((y1-y2)**2))
        positions = np.array((lon_lat1, lon_lat2)) 
        indices2, distances2 = fvcom_pyfvcom.horizontal_transect_nodes(positions)
        for i in range(len(distances2)):
            distances2[i] = true_dist2*distances2[i]/distances2[-1] +true_dist1+true_dist
        print(f"Using outer node, start2={indices2[0]}, end2={indices2[-1]}")
    else:
        lon_lat1, lon_lat2 =(139.795921, 35.2481378) , (139.747317, 35.1262927) 
        x1,y1 =p(lon_lat1[0],lon_lat1[1])[0],p(lon_lat1[0],lon_lat1[1])[1]
        x2,y2 =p(lon_lat2[0],lon_lat2[1])[0],p(lon_lat2[0],lon_lat2[1])[1]
        true_dist2 =math.sqrt((x1-x2)**2 + ((y1-y2)**2))
        positions = np.array((lon_lat1, lon_lat2)) 
        indices2, distances2 = fvcom_pyfvcom.horizontal_transect_nodes(positions)
        for i in range(len(distances2)):
            distances2[i] = true_dist2*distances2[i]/distances2[-1] +true_dist1+true_dist
        print(f"start2={indices2[0]}, end2={indices2[-1]}")


    #print(distances,distances1,distances2)
    indices.extend(indices1);indices.extend(indices2)
    distances_all = np.r_[distances, distances1,distances2]
    #print(distances_all)
    for i,t in enumerate(timestep):
     
        if variable == 'temperature':
            c = fvcom_pyfvcom.data.temp[t, :, indices]
            ## colorbar label
            var = 'Temperature' ; unit = 'degC'  ## Manually
        elif variable == 'salinity':
            c = fvcom_pyfvcom.data.salinity[t, :, indices]
            ## colorbar label
            var = 'salinity' ; unit = 'PSU'  ## Manually
        figsize=(20,9);cmap = 'jet'
        #cmap=cm.balance
        ##### End set parameters

        cb_label = ("{} ({})").format(var, unit)
        png ='./png/vertical/'+var +filename+ '_' +str(timestep_base[i])+ 'profile.png'
        x = distances_all / 1000  # to km from m
        y = fvcom_pyfvcom.grid.siglay_z[:, indices]
       # print(len(x),np.size(y),np.size(c))
        plot = Depth(fvcom_pyfvcom, figsize=figsize, cb_label=cb_label, cmap=cmap)
        ## fill_seabed makes the part of the plot below the seabed gray.
        if variable == 'salinity':
            plot.plot_slice(x, y, c, fill_seabed=True, shading='gouraud',norm=Normalize(vmin=25,vmax=35))
        else:
            plot.plot_slice(x, y, c, fill_seabed=True, shading='gouraud')
        #plot.plot_slice(x, y, c, fill_seabed=True, edgecolors='white')
        plot.axes.set_xlim(right=x.max())  # set the x-axis to the data range
        plot.axes.set_xlabel('Distance (km)')
        plot.axes.set_ylabel('Depth (m)')
        ## Save the figure.
        plot.figure.savefig(png, dpi=300, bbox_inches='tight')
        #print(f"saved figure,{png}")
        plt.close()
        if kh:
            c = fvcom_pyfvcom.data.kh[t, :30, indices]
            var = 'vertical_turbulent_eddy_visicosity' ; unit = 'm^2/s'  ## Manually
            figsize=(20,9);cmap = 'jet'
            cb_label = ("{} ({})").format(var, unit)
            png ='./png/vertical/'+var +filename+ '_' +str(timestep_base[i])+ 'profile.png'
            x = distances_all / 1000  # to km from m
            y = fvcom_pyfvcom.grid.siglay_z[:, indices]
            #print(len(x),np.size(y),np.size(c))
            plot = Depth(fvcom_pyfvcom, figsize=figsize, cb_label=cb_label, cmap=cmap)
            ## fill_seabed makes the part of the plot below the seabed gray.
            plot.plot_slice(x, y, c, fill_seabed=True, shading='gouraud')
            #plot.plot_slice(x, y, c, fill_seabed=True, edgecolors='white')
            plot.axes.set_xlim(right=x.max())  # set the x-axis to the data range
            plot.axes.set_xlabel('Distance (km)')
            plot.axes.set_ylabel('Depth (m)')
            ## Save the figure.
            plot.figure.savefig(png, dpi=300, bbox_inches='tight')
            #print(f"saved figure,{png}")
            plt.close()
        if viscofh:
            c = fvcom_pyfvcom.data.viscofh[t, :30, indices]
            var = 'horizontal_turbulent_eddy_visicosity_or_diffusivity' ; unit = 'm^2/s'  ## Manually
            figsize=(20,9);cmap = 'jet'
            cb_label = ("{} ({})").format(var, unit)
            png ='./png/vertical/'+var +filename+ '_' +str(timestep_base[i])+ 'profile.png'
            x = distances_all / 1000  # to km from m
            y = fvcom_pyfvcom.grid.siglay_z[:, indices]
            #print(len(x),np.size(y),np.size(c))
            plot = Depth(fvcom_pyfvcom, figsize=figsize, cb_label=cb_label, cmap=cmap)
            ## fill_seabed makes the part of the plot below the seabed gray.
            plot.plot_slice(x, y, c, fill_seabed=True, shading='gouraud')
            #plot.plot_slice(x, y, c, fill_seabed=True, edgecolors='white')
            plot.axes.set_xlim(right=x.max())  # set the x-axis to the data range
            plot.axes.set_xlabel('Distance (km)')
            plot.axes.set_ylabel('Depth (m)')
            ## Save the figure.
            plot.figure.savefig(png, dpi=300, bbox_inches='tight')
            #print(f"saved figure,{png}")
            plt.close()
    print('done.')





####################################USER DIFINITION#########################################################################
outer = False
sigmoid = False
kh = True;viscofh=True
####################################END OF USER DIFINITION##################################################################


date_index = pd.date_range("2020-01-01 00:00:00", periods=366, freq="H") #specify start date and output frequency
date_ary = date_index.to_series().dt.strftime("%Y-%m-%d %H:%M:%S")
run = sys.argv[1]
obcfile = f"../{run}/input/input_testcase/TokyoBay_obc.dat" #specify where the obcfile is
variables = [sys.argv[-2],sys.argv[-1]]

if len(sys.argv) <=5:
    filenames = [sys.argv[2]]
else:
    filenames = [sys.argv[i] for i in range(2,len(sys.argv)-2)]

for filename in filenames:

    print(f"begin plotting of {filename}")
    f = '../'+run+'/'+filename+ '_0001.nc'
    fvcom = mod_fvcom.Fvcom2D(ncfile_path=f, obcfile_path=obcfile, m_to_km=True, offset=False)
    stn ={'chiba1buoy':(139.9542517,35.53703833),'chibaharo':(140.0233033,35.61095833),'kawasaki':(139.8340267,35.49019),'urayasu':(139.9417417,35.640085)}
    stn_node = find_closest_node(stn)

    #make_animation(filename)
    
    for variable in variables:
        #if variable ==None:
        #    break
        #for key in stn_node.keys():
        #    index = stn_node[key]
        #    result = plotter2(variable,key,index)
        #    print(result)

        #if outer:
        #  plot_outer_sea(filename)

        plot_vertical(filename)
#plot_map(filenames[0])
























"""
print('start making time contour...')
fname = filename
if sys.argv[2] == 'salinity':
    for place in stn_node:
        index = stn_node[place]
        time = Time(fvcom, figsize=(35,12), cmap = 'jet',cb_label='{} ({})'.format(fvcom.atts.salinity.long_name,
                                                                    fvcom.atts.salinity.units))
        z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
        ## print("z.shape", z.shape)
        ## fill_seabed makes the part of the plot below the seabed grey.
        ## We need to squeeze the data array since we've only extracted a single position.
        ## print(fvcom.grid.h.shape)
        ## print(fvcom.grid.h)
        time.plot_surface(z, fvcom.data.salinity, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h)#,#norm=Normalize(vmin=15, vmax=35))
        time.axes.set_ylabel('Depth (m)')
        time.figure.savefig('./png/contour/salinity'+fname+'_'+place+'.png')
        print('savefig,path=./png/contour/salinity'+fname+'_'+place+'.png')
elif sys.argv[2] == 'temperature':
    for place in stn_node:
        index = stn_node[place]
    #gauge = (140.05, 35.57)  # a sample (lon, lat) position for Estuary example
    #index = fvcom.closest_node(gauge).item()  ### ndarray -> scalar
    # print(index) # 295
        time = Time(fvcom, figsize=(35,12), cmap = 'jet',cb_label='{} ({})'.format(fvcom.atts.temp.long_name,
                                                                    fvcom.atts.temp.units))
        z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
        ## print("z.shape", z.shape)
        ## fill_seabed makes the part of the plot below the seabed grey.
        ## We need to squeeze the data array since we've only extracted a single position.
        ## print(fvcom.grid.h.shape)
        ## print(fvcom.grid.h)
        time.plot_surface(z, fvcom.data.temp, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h)
        time.axes.set_ylabel('Depth (m)')
        time.figure.savefig('./png/contour/temp'+fname+'_'+place+'.png')
        print('savefig,path=./png/contour/temp'+fname+'_'+place+'.png')
"""
######################################################################################################    














    
"""   
##########################################MAKE 3D PLOT HTML##################################################################
    
    
    
# Select item to be plotted.
item='temp'; title = "temprature"; z_range=(15,35)

# Set x, y, and z ranges. = None for auto setting (you may change after plotting).
x_range=None
y_range=None

# Set timestamp_location. = None to deactivate
timestamp_location=(270, 148)

# Create color map.
dmap_ds=fvcom.dmap(item=item, x_range=x_range, y_range=y_range, z_range=z_range, \
                   timestamp_location=timestamp_location)
# Create mesh.
mesh = fvcom.plot_mesh(line_color='black', line_width=0.1)

# Create coastline.
coastline = fvcom.plot_coastline(color='black', line_width=1)

# Overlay as you like
#p = dmap_ds * mesh * coastline

# Apply customization options.
# Customize colorbar.
colorbar_opts=dict(title=title, title_text_font_size="16pt", 
                   title_text_baseline = "bottom",
                   title_standoff =10, title_text_font_style = "normal",
                   major_label_text_font_size = "14pt", major_label_text_align="left",
                   label_standoff=4, width=15)
# Customize Image plot.
p_opts = opts.Image(tools=['hover'], cmap='rainbow', colorbar_position='right',
                    frame_width=400, aspect='equal', colorbar=True, colorbar_opts=colorbar_opts)
#p.redim.label(x='X (km)', y='Y (km)').opts(p_opts)
dmap_ds = dmap_ds.redim.label(x='X (km)', y='Y (km)').opts(p_opts)
# Overlay as you like
p = dmap_ds * mesh * coastline


save = 1
if save == 1:
    dirpath = "./"
    core = "temp_"+f
    outputfile = f"{dirpath}{core}"
    hv.save(p, outputfile, fmt='html', backend='bokeh')

############################make Contour#####################################################################################
#make contor
import numpy as np
## Load an FVCOM model output and plot surface
from PyFVCOM.read import FileReader
from PyFVCOM.read import ncread
from PyFVCOM.plot import Plotter, Time, Depth
from PyFVCOM.tide import make_water_column
from PyFVCOM.grid import unstructured_grid_depths
from cmocean import cm
import matplotlib.pyplot as plt
from pyproj import Proj
import sys
import warnings
# warnings.simplefilter('ignore') # Comment out to show warnings

#fvcom = FileReader(f, dims={'time': slice(0, 1000)}, variables=['zeta', 'temp', 'u'])
#gauge = (140.05, 35.57)  # a sample (lon, lat) position for Estuary example
#index = fvcom.closest_node(gauge).item()  ### ndarray -> scalar
# print(index) # 295

fvcom = FileReader(f, variables=['zeta', 'temp'], dims={'time': slice(0,000), 'node': 150})
time = Time(fvcom, figsize=(40, 9),cmap=cm.thermal, cb_label='{} ({})'.format(fvcom.atts.temp.long_name,
                                                              fvcom.atts.temp.units))
z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
## print("z.shape", z.shape)
## fill_seabed makes the part of the plot below the seabed grey.
## We need to squeeze the data array since we've only extracted a single position.
## print(fvcom.grid.h.shape)
## print(fvcom.grid.h)
time.plot_surface(z, fvcom.data.temp, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h)
time.axes.set_ylabel('Depth (m)')
"""





"""
print('check surface elv at node0...')
fig = plt.figure(figsize=(10,4))
ax = fig.add_subplot(1,1,1,xlabel = 'month',ylabel='salinity')
ax.plot(fvcom.variables['zeta'][:,5],linewidth=0.1,label='observation')
fig.suptitle('sea surface elevation')
ax.legend()
plt.grid()
#plt.ylim(15,35)
plt.xticks([24*30*i for i in range(12)], [i for i in range(1,13)])
#plt.xticks([24*30*i for i in range(12)], [30*i for i in range(12)])
plt.plot()
os.makedirs('./png/zeta',exist_ok=True)
fig.savefig('./png/zeta/'+filename+'_zeta.png')
"""
