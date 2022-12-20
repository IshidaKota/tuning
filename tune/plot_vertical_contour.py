#dev_fvcom_2d_ds.ipynbのデータを元に作成することを想定。鉛直コンターはpyfvcomを用いる。
#pyfvcom
from PyFVCOM.plot import Plotter, Depth
from PyFVCOM.grid import unstructured_grid_depths
from PyFVCOM.read import FileReader
#dev_fvcom_2d_ds
import mod_fvcom
#other
import pandas as pd
from holoviews import opts
import holoviews as hv
import numpy as np
import sys
from pyproj import Proj
import math
#import os


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

print(f"begin plotting of {sys.argv[1]}")
f = '../run/'+sys.argv[1]+ '_0001.nc'
obcfile = '../input_testcase5/TokyoBay_obc.dat' #specify where the obcfile is
variable = sys.argv[2]


fvcom = mod_fvcom.Fvcom2D(ncfile_path=f, obcfile_path=obcfile, m_to_km=True, offset=False)
stn ={'chiba1buoy':(139.9542517,35.53703833),'chibaharo':(140.0233033,35.61095833),'kawasaki':(139.8340267,35.49019),'urayasu':(139.9417417,35.640085)}
stn_node = find_closest_node(stn)

stn_node['miura'] =801
path = './png/map/'
switch = 0
if switch == 1:
    # station map
    node1=stn_node['chiba1buoy']; node2=stn_node['chibaharo'];node3=stn_node['kawasaki'] ;node4=stn_node['urayasu'];node5=stn_node['miura']
    coords1 = [[fvcom.variables['x'][node1].item()/1000, fvcom.variables['y'][node1].item()/1000]]
    coords2 = [[fvcom.variables['x'][node2].item()/1000, fvcom.variables['y'][node2].item()/1000]]
    coords3 = [[fvcom.variables['x'][node3].item()/1000, fvcom.variables['y'][node3].item()/1000]]
    coords4 = [[fvcom.variables['x'][node4].item()/1000, fvcom.variables['y'][node4].item()/1000]]
    coords5 = [[fvcom.variables['x'][node5].item()/1000, fvcom.variables['y'][node5].item()/1000]]
    node6 =46;node7=531;node8=637,;node9=930
    l1 =(fvcom.variables['x'][node6].item()/1000, fvcom.variables['y'][node6].item()/1000)
    l2 =(fvcom.variables['x'][node7].item()/1000, fvcom.variables['y'][node7].item()/1000)
    l3 =(fvcom.variables['x'][node8].item()/1000, fvcom.variables['y'][node8].item()/1000)
    l4 =(fvcom.variables['x'][node9].item()/1000, fvcom.variables['y'][node9].item()/1000)
    #opts.defaults(opts.Points(tools=['hover']))
    path1 = hv.Path([l1,l2])
    path2 = hv.Path([l2,l3])
    path3 = hv.Path([l3,l4])
    mesh = fvcom.plot_mesh(line_color='blue', line_width=0.1)
    points1 = fvcom.plot_point(coords1, color='#63000f', marker='x', size=12,label='chiba1buoy')
    points2 = fvcom.plot_point(coords2, color='green', marker='x', size=12,label='chibaharo')
    points3 = fvcom.plot_point(coords3, color='red', marker='x', size=12,label='kawasaki')
    points4 = fvcom.plot_point(coords4, color='blue', marker='x', size=12,label='urayasu')
    points5 = fvcom.plot_point(coords5, color='blue', marker='x', size=12,label='miura')
    p = mesh * points1 * points2 * points3 * points4 * points5 * fvcom.plot_coastline(color='black',linewidth=1) *path1 *path2 * path3
    p_opts = opts.Points(frame_width=500, aspect='equal', xlim=(340, 460))
    p.opts(p_opts)

    
    hv.save(p, path+'st_location.png', fmt='png',dpi=300, backend='bokeh')

##########################
"""
平面図
# Set PNG output file path
#dirpath = "./png/horizontal/"
#core = "tri_temp"+sys.argv[1]
# toolbar='right' for interactive view in a browser
#     Possible to zoom in on the browser and output png in it.
# toolbar=None to deactivate this interactivity. Useful for PDF export using pyppeteer.
#toolbar=None # None, 'right'
toolbar = 'right'
# Select item to be plotted.
if sys.argv[2] == 'salinity':
    item='salinity'; title = "Salinity (psu)"; z_range=(15,35)
else:
    item='temp'; title = "temperature"; z_range=(0,30)

# Set time and sigma level
timestep=[210,240,270]
sigmalay=[2,15,28]

# Set x and y ranges. = None for auto setting (you may change after plotting).
x_range=None#(340,450) # None, (min, max)
y_range=None#(3820,3960) # None, (min, max)

# Set timestamp_location. = None to deactivate
timestamp_location=(380, 3920)



colorbar_opts=dict(title=title, title_text_font_size="16pt", 
                   title_text_baseline = "bottom",
                   title_standoff =10, title_text_font_style = "normal",
                   major_label_text_font_size = "14pt", major_label_text_align="left",
                   label_standoff=4, width=15)
p_opts = opts.Image(tools=['hover'],cmap='rainbow', colorbar_position='right',
                    frame_width=500, aspect='equal', colorbar=True, colorbar_opts=colorbar_opts,
                    toolbar=toolbar)
                    #hooks=[fixBottomMargin])

for time in timestep:
    for sigma in sigmalay:
        #outputfile = f"{dirpath}{core}{time:03}{sigma:03}"

        # Create color map.
        smap = fvcom.splot_2dh(item=item, time=time, sigma=sigma,
                            x_range=x_range, y_range=y_range, z_range=z_range,
                            timestamp_location=timestamp_location)
        smap=smap.redim.label(x='X (km)', y='Y (km)').opts(p_opts)

        # Create mesh.
        mesh = fvcom.plot_mesh(line_color='black', x_range=x_range, y_range=y_range,line_width=0.2)

        # Create coastline.
        coastline = fvcom.plot_coastline(color='black',x_range=x_range, y_range=y_range, line_width=1)

        # Overlay as you like
        p = smap * mesh * coastline
        hv.save(p,'./png/horizontal/'+item+sys.argv[1]+'t'+str(time)+'siglay'+str(sigma)+'.png', fmt='html',backend='bokeh')


"""




#dt = datetime.datetime.fromtimestamp(p.stat().st_ctime)

########################plot vertical contour####################################################
fvcom = FileReader(f, variables=['zeta', 'salinity','temp'], dims={'time': slice(0,8700)})
##### Start set parameters
timestep = [30,90,150,210,240,270]
#timestep = [240]
    # lon_lat1, lon_lat2 = (138, 33), (139, 33)     ## Miwa's mesh
lon_lat1, lon_lat2 = (140.01255, 35.6169609), (139.71812, 35.3276445) ## Estuary example
positions = np.array((lon_lat1, lon_lat2)) 
indices, distances = fvcom.horizontal_transect_nodes(positions)

lon_lat1, lon_lat2 = (139.71812, 35.3276445), (139.800387, 35.2118708) ## Estuary example
positions = np.array((lon_lat1, lon_lat2)) 
indices1, distances1 = fvcom.horizontal_transect_nodes(positions)
for i in range(len(distances1)):
    distances1[i] += distances[-1]

lon_lat1, lon_lat2 = (139.800387, 35.2118708), (139.561148, 34.9360547) ## Estuary example
positions = np.array((lon_lat1, lon_lat2)) 
indices2, distances2 = fvcom.horizontal_transect_nodes(positions)
for i in range(len(distances2)):
    distances2[i] += distances1[-1]

indices.extend(indices1);indices.extend(indices2)
distances = np.r_[distances, distances1,distances2]
for t in timestep:

    if variable == 'temperature':
        c = fvcom.data.temp[t, :, indices]
        ## colorbar label
        var = 'Temperature' ; unit = 'degC'  ## Manually
    else:
        c = fvcom.data.salinity[t, :, indices]
        ## colorbar label
        var = 'salinity' ; unit = 'PSU'  ## Manually
    figsize=(20,9);cmap = 'jet'
    #cmap=cm.balance
    ##### End set parameters

    cb_label = ("{} ({})").format(var, unit)
    png ='./png/vertical/'+var +sys.argv[1]+ '_' +str(t)+ 'profile.png'
    x = distances / 1000  # to km from m
    y = fvcom.grid.siglay_z[:, indices]
    #print(type(x),type(y),x,y)
    plot = Depth(fvcom, figsize=figsize, cb_label=cb_label, cmap=cmap)
    ## fill_seabed makes the part of the plot below the seabed gray.
    plot.plot_slice(x, y, c, fill_seabed=True, shading='gouraud')
    #plot.plot_slice(x, y, c, fill_seabed=True, edgecolors='white')
    plot.axes.set_xlim(right=x.max())  # set the x-axis to the data range
    plot.axes.set_xlabel('Distance (km)')
    plot.axes.set_ylabel('Depth (m)')
    ## Save the figure.
    plot.figure.savefig(png, dpi=600, bbox_inches='tight')