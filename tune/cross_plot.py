#plot different simulations at the same time.
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

places = ['chiba1buoy','chibaharo','kawasaki','urayasu']
variables = ['salinity','temperature']
units = ['PSU','degC']

for variable,unit in zip(variables,units):
    print(variable,unit)
    for place in places:
        fig=plt.figure(figsize=(6,3))
        ax=fig.add_subplot(1,1,1,ylabel=variable+'('+unit+')',xlabel='month(2020)')
        paths =glob('./df/result/'+place+variable+'*')
        for path in paths:
            df = pd.read_csv(path)
            #print(df.head())
            casename = path[-10:-4]
            print(casename)
            ax.plot(df['2'][:],label='2_'+casename)
            ax.plot(df['15'][:],label='15_'+casename)
            ax.plot(df['28'][:],label='28_'+casename)
        ax.legend();ax.grid();plt.plot()
        plt.suptitle('cross plotting of '+variable+' of '+place)
        plt.xticks([0,30,59,90,120,150,181,212,243,273,304,334],[i for i in range(1,13)])
        fig.savefig('./png/cross_plot/'+variable+'_'+place+'.png',\
            bbox_inches='tight', pad_inches=0.1,dpi=600)


