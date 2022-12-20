import netCDF4
import sys
import matplotlib.pyplot as plt
import datetime
from jdcal import  jd2gcal ,MJD_0
from matplotlib import dates as mdates
Filenames = [sys.argv[i] for i in range(1,len(sys.argv))]
print(Filenames)
Ncs = [netCDF4.Dataset('../run27/'+filename+ '_0001.nc','r') for filename in Filenames ]
dates_float = [jd2gcal(MJD_0,Ncs[0].variables['Itime'][i]) for i in range(len(Ncs[0].variables['Itime']))] #ymd
dates = [datetime.datetime(dates_float[i][0],dates_float[i][1],dates_float[i][2],\
        hour=int(Ncs[0].variables['Itime2'][i]/3600000))for i in range(len(Ncs[0].variables['Itime']))] #ymd+h
print(Ncs[0].variables.keys())
Variables = ['viscofh','kh']
HPRNU = ['0.01','0.1','1.0','10.0']
index = 920

def plot_viscof(Variables,sigma):
    if sigma ==2:
        layer = 'upper'
    if sigma ==28:
        layer = 'bottom'
    fig,axs = plt.subplots(len(Variables),1,figsize=(10,10))
    for i,nc in enumerate(Ncs):
        for j,variable in enumerate(Variables):
            if variable in nc.variables.keys():
                val_upper = nc.variables[variable][:,sigma,index]
             
                if variable =='kh' or variable == 'viscofh':
                    val_upper = val_upper*(10**4)
                    #print(type(dates),type(val_upper))
                    axs[j].plot(dates,val_upper,label=layer+'_'+HPRNU[i])
                    axs[j].set_ylabel(variable + ' (cm^2/s)')
                else:
                    axs[j].plot(dates,val_upper,label=layer+'_'+ HPRNU[i])
                    axs[j].set_ylabel(variable + ' (m^2/s)')
                axs[j].xaxis.set_major_formatter(mdates.DateFormatter('%y-%m'))
                axs[j].legend()
                axs[j].set_title(variable+' '+layer+ ' at chiba1buoy')

                if variable == 'kh':
                    axs[j].set_ylim(0,5)
                if variable == 'viscofh':
                    axs[j].set_ylim(0,4*10**6)
    plt.tight_layout()
    fig.savefig('./png/viscof/'+Variables[0]+layer+Filenames[0]+'to'+Filenames[-1]+'.png'\
                ,bbox_inches='tight', pad_inches=0.1,dpi=300)
    plt.close()
plot_viscof(['kh','viscofh'],2)
plot_viscof(['km','viscofm'],2)
plot_viscof(['kh','viscofh'],28)
plot_viscof(['km','viscofm'],28)