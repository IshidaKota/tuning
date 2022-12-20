#plot different simulations at the same time.
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import sys
import matplotlib.colors as mcolors
import datetime
from jdcal import gcal2jd, jd2gcal ,MJD_0
from matplotlib import dates as mdates
places = ['chiba1buoy','chibaharo','kawasaki','urayasu']
variables = ['salinity','temperature']
units = ['PSU','degC']
files = [sys.argv[i] for i in range(1,len(sys.argv))]
num = "".join([file[-1] for file in files])
#cmaps = plt.get_cmap("tab10")
#print(cmaps)
css_colors = mcolors.BASE_COLORS
linestyles = ['solid','dotted','dashed','dashdot']

print(files)
for variable,unit in zip(variables,units):
    print(variable,unit)
    for place in places:
        fig=plt.figure(figsize=(6,4))
        ax=fig.add_subplot(1,1,1,ylabel=variable+' ('+unit+')',xlabel='month(2020)')
        
        for file,color in zip(files,css_colors):
            path ='./df/result/'+place+variable+'_'+file+'_.csv'
            df = pd.read_csv(path)
            column = pd.to_datetime(df['date'])
    
            ax.plot(df['2'][:],label='2_'+file,color=color,linestyle=linestyles[0])
            #ax.plot(df['15'][:],label='15_'+file,color=color,linestyle=linestyle)
            ax.plot(df['28'][:],label='28_'+file,color=color,linestyle=linestyles[1])
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))
        ax.legend();ax.grid();plt.plot()
        plt.suptitle('cross plotting of '+variable+' of '+place)
        plt.xticks([0,30,59,90,120,150,181,212,243,273,304,334],[i for i in range(1,13)])
        fig.savefig('./png/cross_plot/'+variable+'_'+place+num+'.png',\
            bbox_inches='tight', pad_inches=0.1,dpi=300)
