import sys
import os
from functools import partial, partialmethod
import numpy as np
import pandas as pd
import pyproj
from pyproj import Proj
### import xarray as xr
import netCDF4
### import matplotlib.pyplot as plt
import matplotlib.tri as tri
import holoviews as hv
#import datashader as ds
import datashader.utils as du
#import datashader.transfer_functions as tf
from holoviews.operation.datashader import rasterize

#from holoext.xbokeh import Mod
import geoviews as gv

import subprocess

#import pdfkit
#from holoviews.plotting.bokeh.element import (line_properties, fill_properties, text_properties)
hv.extension('bokeh', 'matplotlib')

def proj_trans(df=None, col_x='Longitude', col_y='Latitude', epsg0="2451", epsg1="4612"):
    '''
    Transform coordinates of pandas.DataFrame (first two columns of (x, y)) 
    from epsg0 to epsg1; then the column names of (x, y) are renamed to (col_x, col_y),
    and return the DataFrame
    
    Parameters
    -------------
    df : pandas.DataFrame
    col_x : str
        x coordinate name
    col_y : str
        y coordinate name
    epsg0 : str
        original epsg code
    epsg1 : str
        transformed epsg code

    returns
    -------
    transformed pandas.DataFrame
    '''

    EPSG0 = pyproj.Proj("+init=EPSG:" + epsg0)
    EPSG1 = pyproj.Proj("+init=EPSG:" + epsg1)
    #x1,y1 = pyproj.transform(EPSG0, EPSG1, df.iloc[:,0].values, df.iloc[:,1].values)
    p = Proj(proj='utm',zone=54,ellps='WGS84', preserve_units=False)
    x1,y1 = p(df.iloc[:,0].values,df.iloc[:,1].values)
    ### print(len(x1))
    ### print(len(df))
    df[df.columns[0]]=x1
    df[df.columns[1]]=y1
    df.rename(columns = {df.columns[0]:col_x, df.columns[1]:col_y}, inplace=True)
    return df

class Fvcom2D:
    '''
    Postprocessing for FVCOM 4 netcdf output in its original coordinates
    '''
    def __init__(self, ncfile_path='tst_0001.nc', 
                 obcfile_path='tst_obc.dat', m_to_km=True, offset=False):
        '''
        Parameters
        ----------
        ncfile_path : str
            FVCOM output netcdf file
        obcfile_path : str
            FVCOM casename_obc.dat file
        m_to_km : bool
            = True if converting x and y axis units from m to km.
        offset : bool
           offset = True if the left-bottom corner is set to the origin.

        Returns
        -------
        Instance of class FVCOM2D
        '''

        self.fvcom_nc=ncfile_path
        if os.path.isfile(self.fvcom_nc):
            self.FVCOM = netCDF4.Dataset(self.fvcom_nc, 'r')
        else:
            print(f"ERROR: File {self.fvcom_nc} does not exit.")
            sys.exit()
        self.variables = self.FVCOM.variables
        #print(f"Dictionaly keys of netCDF4.Dataset = {self.variables.keys()}")
        self.fobc=obcfile_path # = None if not exist
        if os.path.isfile(self.fobc):
            df = pd.read_csv(self.fobc, header=None, skiprows=1, delim_whitespace=True)
            ### -1 because index in FVCOM starts from 1 while from 0 in Python
            self.node_bc = df.iloc[:,1].values - 1
            #print(f"Open boundary nodes = {self.node_bc}")
        else:
            print(f"ERROR: File {self.fobc} does not exit.")
            sys.exit()            
        self.x, self.y = self.variables['x'][:], self.variables['y'][:] ### at node
        self.xc, self.yc = self.variables['xc'][:], self.variables['yc'][:]   ### at cell center
        # Convert axis units from m to km
        if m_to_km:
            self.x /= 1000.0; self.y /= 1000.0; self.xc /= 1000.0; self.yc /= 1000.0
        # Set the left-bottom corner as the coordinate origin. 
        if offset:
            xmin = self.x.min(); ymin = self.y.min()
            self.x -= xmin; self.y -= ymin; self.xc -= xmin; self.yc -= ymin
        # Gets sigma coordinate values
        self.siglay, self.siglev = self.variables['siglay'][:], self.variables['siglev'][:]
        # Get time variables
        self.iint = self.variables['iint'][:]
        self.time = self.variables['time'][:]
        self.Itime = self.variables['Itime'][:]
        self.Itime2 = self.variables['Itime2'][:]
        # Gets bathymetry
        self.z = self.variables['h'][:]
        # Gets verts = [(x0,y0,h0),[x1,y1,h1],...,[xn,yn,hn]]
        verts = [(xi, yi, hi) for xi, yi, hi in zip(self.x.data, self.y.data, self.z.data)]
        self.verts = pd.DataFrame(verts, columns=['x', 'y', 'z'])
        # Gets connectivity array nv (3 node numbers of element i)
        nv = self.variables['nv'][:].T - 1
        self.triang = tri.Triangulation(self.x, self.y, triangles=nv)
        # Since node indexes in an element are defined clockwise in FVCOM,
        #     change them to counterclockwise, which is the matplotlib.tri specification.
        # FVCOMでは時計回りに定義されているので，matplotlib.triの仕様である反時計回りに変更する．
        # node番号自体は不変．
        self.nv=nv[:,::-1]
        self.tris = pd.DataFrame(self.nv, columns=['v0', 'v1', 'v2'])
        self.mesh = du.mesh(self.verts, self.tris) 
        self.trimesh = hv.TriMesh((self.tris, self.verts))
        # Element index of 3 neighbors of element i
        # 各セルには一般に3つのセルが隣接（境界セルはこの限りではない）
        self.nbe = np.array([[self.nv[n, j], self.nv[n, (j+2)%3]] \
            for n in range(len(self.triang.neighbors)) \
            for j in range(3) if self.triang.neighbors[n,j] == -1])

    @property
    def timestamp(self):
        '''
        Create time stamp list in ["D HH:MM:SS"]
        '''
        dayi = self.Itime
        hourf = self.Itime2 / 3600000
        houri = hourf.astype('int32')
        minf = (hourf - houri) * 60
        mini = minf.astype('int32')
        secf = (minf - mini) * 60
        seci = secf.astype('int32')
        return [f"{d} {h:02}:{m:02}:{s:02}" for d, h, m, s in zip(dayi, houri, mini, seci)]

    def plot_mesh(self, x_range=None, y_range=None, **kwargs):
        '''
        Plot mesh
        
        Parameters
        ----------
        x_range : tuple
            Plot range of x axis (xmin, xmax)
        y_range : tuple
            Plot range of y axis (ymin, ymax)
        '''

        # Gets **kwargs with default values
        #width=kwargs.get('width', 300)
        #height=kwargs.get('height', 300)
        line_color=kwargs.get('line_color', 'blue')
        line_width=kwargs.get('line_width', 0.1)

        # jsasaki debug
        #p = self.trimesh.edgepaths.options(width=width, height=height, 
        #                                   line_width=line_width, line_color=line_color)
        p = self.trimesh.edgepaths.options(line_width=line_width, line_color=line_color)

        # Set plot range
        if x_range is None and y_range is None:
            return p
        elif y_range is None:
            return p.redim.range(x=x_range)
        elif x_range is None:
            return p.redim.range(y=y_range)
        else:
            return p.redim.range(x=x_range, y=y_range)
    
    def plot_coastline(self, x_range=None, y_range=None, **kwargs):
        '''
        Plot coastline. May take much time for interactive plotting.

        Parameters
        ----------
        x_range : tuple
            Plot range in x-axis (xmin, xmax)
        y_range : tuple
            Plot range in y-aixs (ymin, ymax)
        '''

        # Gets **kwargs with default values
        #width=kwargs.get('width', 300)
        #height=kwargs.get('height', 300)
        color=kwargs.get('color', 'red')
        line_width=kwargs.get('line_width', 0.5)

        # paths = [[(x0s, y0s), (x0e, y0e)], ..., [(xns, yns), (xne, yne)]]
        paths = [[(self.x[self.nbe[m,:]][0], self.y[self.nbe[m,:]][0]), \
                  (self.x[self.nbe[m,:]][1], self.y[self.nbe[m,:]][1])] \
                 for m in range(len(self.nbe))]

        # jsasaki debug
        #p = hv.Path(paths).options(width=width, height=height, 
        #                           color=color, line_width=line_width)
        p = hv.Path(paths).options(color=color, line_width=line_width)

        if x_range is None and y_range is None:
            return p
        elif y_range is None:
            return p.redim.range(x=x_range)
        elif x_range is None:
            return p.redim.range(y=y_range)
        else:
            return p.redim.range(x=x_range, y=y_range)

    def get_val(self, item, time=None, sigma=None):
        '''
        Gets item values dependent considering its shape with specifying
        output time step (if time sereis data) and sigma level.

        Parameters
        ----------
        item : str
            Name of item in dictionary keys
        time : int
            Time series index
        sigma : int
            Layer number in sigma coordinates

        Returns
        -------
        pandas.DataFrame
        '''

        # Check whether item exists among the keys
        if not item in self.variables.keys():
            print(f"Error: Item {item} does not exit in keys.")
            sys.exit()
        val = self.FVCOM.variables[item]

        # Check the item dimensions
        if len(val.shape) == 1 : # 1-D item
            scalar = val
        elif len(val.shape) >= 2: # 2-D or 3-D item
            if time < 0 or time > val.shape[0]:
                print(f"ERROR: Time index = {time} is out of range.")
                print(f"Time index should be between 0 and {val.shape[0]-1}.")
                sys.exit()
            if len(val.shape) == 2: ### 2-D item
                scalar = val[time]
            elif len(val.shape) == 3: ### 3-D item
                if sigma < 0 or sigma > val.shape[1]:
                    print(f"ERROR: Sigma layer number = {sigma} is out of range.")
                    print(f"Sigma layer number should be between 0 and {val.shape[1]-1}.")
                    sys.exit()
                else:
                    scalar = val[time][sigma]
            else:
                print("ERROR: the shape is incorrect.")
                sys.exit()

        verts = [(xi, yi, hi) for xi, yi, hi in zip(self.x.data, self.y.data, scalar)]
        return pd.DataFrame(verts, columns=['x', 'y', item])

    def __plot_2dh(self, item, time, sigma):
        '''
        Internal method for converting pd.DataFrame of item to hv.TriMesh

        Parameters
        ----------
        item : str
        time : int
            Time index to be read
        sigma :int
            Sigma level index to be read
        Returns
        -------
        HoloViews.TriMesh
        '''

        # Reading pd.DataFrame of item using get.val() method 
        if len(self.variables[item].shape) == 1:   ### 1-D item (e.g.: depth) 
            df_item = self.get_val(item)
        elif len(self.variables[item].shape) == 2: ### 2-D item (e.g.: zeta)
            df_item = self.get_val(item, time=time)
        elif len(self.variables[item].shape) == 3: ### 3-D item (e.g.: salainity)
            df_item = self.get_val(item, time=time, sigma=sigma)
        else:
            print('ERROR: The shape of item is incorrect in self.__plot_2dh')
            sys.exit()
        return hv.TriMesh((self.tris, hv.Points(df_item)))

    def dmap(self, item, x_range=None, y_range=None, z_range=None, timestamp_location=None):
        '''
        Interactive plotting using DynamicMap
        Internal functions need to be deifned in order to be an argument of 
            functools.partial().
        Tips: Prepare an internal function,
            and set time and sigma at the end of the arguments in partial.
        内部関数を用意すること，partialの引数の最後をtimeとsigmaにするのがポイント
        
        Parameters
        ----------
        self : FVCOM2D instance
        item : str
        x_range : tuple
        y_range : tuple, optional
        z_range : tuple, optional
        timestamp_location : tuple, optional
        '''

        def item3d_plot(item, timestamp_location, sigma, time):
            '''
            Internal function for plotting 3D item using functools.partial().
            partialの引数とするため，内部関数とする.

            Parameters
            ----------
            item : str
            timestamp_location : tuple
            sigma : int
                sigma level
            time : int
                output timestep of time series data

            Returns
            -------
            hv.TriMesh * self.plot_timestamp
            '''
            #item1 = self.get_val('u',time=time,sigma=sigma)
            #item2 = self.get_val('v',time=time,sigma=sigma)
            #x,y  = np.mgrid[-10:10,-10:10] * 0.25
            #sine_rings  = np.sin(item1**2+item2**2)*np.pi+np.pi
            #exp_falloff = item1**2 + item2**2

            #vector_data = (item1,item2, sine_rings, exp_falloff)
            #timestamp_txt = self.plot_timestamp(time, timestamp_location)
            #return hv.TriMesh((self.tris, hv.VectorField(vector_data).opts(magnitude='Magnitude'))) * timestamp_txt

            item = self.get_val(item, time=time, sigma=sigma)
            timestamp_txt = self.plot_timestamp(time, timestamp_location)
            return hv.TriMesh((self.tris, hv.Points(item))) * timestamp_txt

        def item2d_plot(item, timestamp_location, time):
            '''
            Internal function for plotting 2D item (e.g. zeta) using functools.partial().
            
            Parameters
            ----------
            item : str
            timestamp_location : tuple
            time : int
            
            Returns
            -------
            hv.TriMesh * self.plot_timestamp
            '''

            item = self.get_val(item, time=time)
            timestamp_txt = self.plot_timestamp(time, timestamp_location)
            return hv.TriMesh((self.tris, hv.Points(item))) * timestamp_txt        
        
        if len(self.variables[item].shape) == 1:   ### 1-D item (ex: depth) 
            print(f"ERROR: 1-D item of {item} is not supported in FVCOM2D.dmap")
            sys.exit()
        elif len(self.variables[item].shape) == 2: ### 2-D item (ex: zeta)
            p = partial(item2d_plot, item, timestamp_location)
            dmap = hv.DynamicMap(p, kdims=['time']).redim.values(time=range(len(self.time))) #Dynamicの引数を合わせる
        elif len(self.variables[item].shape) == 3: ### 3-D item (ex: salainity)
            p = partial(item3d_plot, item, timestamp_location)
            dmap = hv.DynamicMap(p, kdims=['sigma', 'time']).redim.values( \
                sigma=range(len(self.siglay))).redim.values(time=range(len(self.time)))
        else:
            print('ERROR: The shape of the item is incorrect in FVCOM2D.dmap')
            sys.exit()

        if x_range is None:
            x_range = (self.x.min(), self.x.max())
        if y_range is None:
            y_range = (self.y.min(), self.y.max())
        if z_range is None:
            return rasterize(dmap, x_range=x_range, y_range=y_range)
        else:
            return rasterize(dmap.redim.range(**{item : z_range}), \
                             x_range=x_range, y_range=y_range)

    def splot_2dh(self, item, time=0, sigma=0, x_range=None, y_range=None,
                  z_range=None, timestamp_location=None, **kwargs):
        '''
        Static horizontal 2D plot (one sigma and time frame)
        静的な平面2次元プロット．InteractiveにはDynamicMapを用いたdmap()を使用する
        '''
        if x_range is None:
            x_range = (self.x.min(), self.x.max())
        if y_range is None:
            y_range = (self.y.min(), self.y.max())

        p = self.__plot_2dh(item=item, time=time, sigma=sigma)
        if z_range is None or z_range == False:
            p = rasterize(p, x_range=x_range, y_range=y_range)
        else:
            p = rasterize(p.redim.range(**{item : z_range}), x_range=x_range, y_range=y_range)

        timestamp_txt = self.plot_timestamp(time, timestamp_location)
        return p * timestamp_txt

    def plot_point(self, coords, color='red', marker='+', size=10, **kwargs):
        p = hv.Points(coords, **kwargs)
        
        return p.opts(color=color, marker=marker, size=size)
        
    def plot_timestamp(self, time=0, timestamp_location=None, **kwargs):
        '''
        timestamp_location: tuple (x, y)
        Currently setting the location using graphic space unsupported in Holoviews
        スクリーン座標系で位置が設定できるとよいが，現在はHoloviewsの仕様で不可
        timestamp_location=Noneのときは空文字を出力することで汎用化した
        '''
        if timestamp_location is None:
            return hv.Text(0, 0, "", **kwargs)
        else:
            timestamp = self.timestamp[time]
            return hv.Text(timestamp_location[0], timestamp_location[1], timestamp, **kwargs)

    def make_frame(self, item, title=None, time=None, sigma=None,
                   plot_mesh=True, plot_coastline=True,
                   x_range=None, y_range=None, z_range=None,
                   timestamp_location=None):
        '''
        Make a frame at one time step [and one sigma layer (only for 3D)]
        1つのframe出力．アニメーション作成に有効
        '''
        ptri = self.splot_2dh(item=item,time=time, sigma=sigma, colorbar_title=title,
                              x_range=x_range, y_range=y_range,
                              z_range=z_range, timestamp_location=timestamp_location)
        pmesh = self.plot_mesh(x_range=x_range, y_range=y_range)
        pcoastline = self.plot_coastline(x_range=x_range, y_range=y_range)
        #timestamp_txt = self.plot_timestamp(time, timestamp_location)

        if plot_mesh and plot_coastline:
            return ptri * pmesh * pcoastline
        elif plot_mesh:
            return ptri * pmesh
        elif plot_coastline:
            return ptri * pcoastline
        else:
            return ptri

class Fvcom2D_trans(Fvcom2D):
    '''
    Postprocessing for FVCOM 4 with converting coordinates
    Requires function proj_trans()
    (x, y)から(lon, lat)への変換を想定しているが，epsg0とepsg1と座標名を与えることで任意の変換が可能
    '''
    def __init__(self, col_x='Longitude', col_y='Latitude', epsg0="2451", epsg1="4612", \
                 ncfile_path='./output/tokyo_0001.nc', obcfile_path='./inputfile/tokyo_obc.dat'):
        '''
        m_to_km and offset must be False in super().__init__().
        '''
        super().__init__( ncfile_path, obcfile_path, m_to_km=False, offset=False)
        self.verts = proj_trans(df=self.verts, col_x=col_x, col_y=col_y, epsg0=epsg0, epsg1=epsg1)
        points = gv.operation.project_points(gv.Points(self.verts, vdims=['z']))
        self.mesh = du.mesh(self.verts, self.tris)
        self.trimesh = hv.TriMesh((self.tris, points))
        self.verts = points
