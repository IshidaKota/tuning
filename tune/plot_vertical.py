

import numpy as np
import math


from pyproj import Proj
### import xarray as xr


#hv.extension('bokeh', 'matplotlib')


from PyFVCOM.plot import Depth
from PyFVCOM.grid import Domain,find_connected_nodes
from PyFVCOM.read import FileReader
from PyFVCOM.current import vector2scalar


def plot_vertical():
    print('start making vertical contour...')
    f = '../run27/Tokyo5_0001.nc'
    t = 5
    fvcom_pyfvcom = FileReader(f, variables=['zeta', 'salinity','temp','u','v'], dims={'time': slice(0,10)})
    ##### Start set parameters
    #timestep_base = [30,90,150,210,240,270];scale_factor = fvcom.time[1]-fvcom.time[0]
    #timestep = [int(timestep_base[i]*scale_factor) for i in range(len(timestep_base))]
    
    #outer 45 484 589 987
    mesh = Domain(f,'cartesian','54')
    lon_lat1, lon_lat2 = (140.026631, 35.6347797), (140.0, 35.65) 
    positions = np.array((lon_lat1, lon_lat2)) 
    indices, distances = fvcom_pyfvcom.horizontal_transect_elements(positions)
    nodes = [mesh.grid.triangles[i] for i in indices]

    #補正
    p = Proj(proj='utm',zone=54,ellps='WGS84', preserve_units=False)
    x1,y1 =p(lon_lat1[0],lon_lat1[1])[0],p(lon_lat1[0],lon_lat1[1])[1]
    x2,y2 =p(lon_lat2[0],lon_lat2[1])[0],p(lon_lat2[0],lon_lat2[1])[1]
    true_dist =math.sqrt((x1-x2)**2 + ((y1-y2)**2))
    for i in range(len(distances)):
        distances[i] = true_dist*distances[i]/distances[-1]

    

    c ,dist=vector2scalar(fvcom_pyfvcom.data.u[t,:,indices],fvcom_pyfvcom.data.v[t,:,indices])
    
        ## colorbar label
    var = 'velocity_magnitude' ; unit = 'm/s'  ## Manually
    figsize=(20,9);cmap = 'jet'
    #cmap=cm.balance
    ##### End set parameters

    cb_label = ("{} ({})").format(var, unit)
    png ='./png/vertical/velocity_profile.png'
    x = distances / 1000  # to km from m
    Y = []
    for i in range(len(nodes)):
        y=0
        for j in range(len(nodes[i])):
            y += fvcom_pyfvcom.grid.siglay_z[:, nodes[i][j]]
        Y.append(y/len(nodes[i]))
    Y = np.array(Y)
    print(Y,c)
    plot = Depth(fvcom_pyfvcom, figsize=figsize, cb_label=cb_label, cmap=cmap)
    ## fill_seabed makes the part of the plot below the seabed gray.
    plot.plot_slice(x, Y, c, fill_seabed=True, shading='ground')
    #plot.plot_slice(x, y, c, fill_seabed=True, edgecolors='white')
    plot.axes.set_xlim(right=x.max())  # set the x-axis to the data range
    plot.axes.set_xlabel('Distance (km)')
    plot.axes.set_ylabel('Depth (m)')
    ## Save the figure.
    plot.figure.savefig(png, dpi=600, bbox_inches='tight')
    #print(f"saved figure,{png}")
    print('done.')


plot_vertical()