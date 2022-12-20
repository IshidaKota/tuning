## Load an FVCOM model output and plot surface
from PyFVCOM.read import FileReader
from PyFVCOM.plot import Time
from PyFVCOM.grid import unstructured_grid_depths
import cmocean as cm 
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pyproj import Proj
import sys
import math
import mod_fvcom

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
        print(f"Node {res} is {min_dist}km to the station:{key}")
    return index



f = '../run27/'+sys.argv[1]+ '_0001.nc'
obcfile = '../run27/input/input_testcase/TokyoBay_obc.dat' 
fvcom = mod_fvcom.Fvcom2D(ncfile_path=f, obcfile_path=obcfile, m_to_km=True, offset=False)

stn ={'chiba1buoy':(139.9542517,35.53703833),'chibaharo':(140.0233033,35.61095833),'kawasaki':(139.8340267,35.49019),'urayasu':(139.9417417,35.640085)}
stn_node = find_closest_node(stn)
stn = {'chiba1buoy':133,'chibaharo':42,'kawasaki':205,'urayasu':40,'check':366}

##additional
#stn_node = {}
#stn_node['low_temperature'] =  1282

fname = sys.argv[1]
ncfile = "../run27/"+fname + "_0001.nc"


for place in stn_node:
    index = stn_node[place]
    fvcom = FileReader(ncfile, variables=['zeta', 'salinity'], dims={'time': slice(0,8700), 'node': index})
    time = Time(fvcom, figsize=(35,12), cmap = 'jet',cb_label='{} ({})'.format(fvcom.atts.salinity.long_name,
                                                                fvcom.atts.salinity.units))
    z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
#dynamic
    time.plot_surface(z, fvcom.data.salinity, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h,norm=Normalize(vmin=25, vmax=34))
    time.axes.set_ylabel('Depth (m)')
    time.figure.savefig('./png/contour/salinitydyn'+fname+'_'+place+'.png')



    fvcom = FileReader(ncfile, variables=['zeta', 'temp'], dims={'time': slice(0,8700), 'node': index})
    time = Time(fvcom, figsize=(35,12), cmap = 'jet',cb_label='{} ({})'.format(fvcom.atts.temp.long_name,
                                                                fvcom.atts.temp.units))
    z = unstructured_grid_depths(fvcom.grid.h, fvcom.data.zeta, fvcom.grid.siglay)
    time.plot_surface(z, fvcom.data.temp, fill_seabed=True, h_offset=0, h_min=fvcom.grid.h)
    time.axes.set_ylabel('Depth (m)')
    time.figure.savefig('./png/contour/tempdyn'+fname+'_'+place+'.png')
    print('savefig,path=./png/contour/tempdyn'+fname+'_'+place+'.png')