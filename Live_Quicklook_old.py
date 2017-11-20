# -*- coding: utf-8 -*-
"""
Program to produce Quicklooks from the MBR2 Cloud Radar and 2m meteorology data
Marcus Klingebiel, Tobias Machnitzki, April 2017
(marcus.klingebiel@mpimet.mpg.de)
(tobias.machnitzki@mpimet.mpg.de)
"""
#%%
#=========================
# # Import libraries:
#=========================

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
from matplotlib import gridspec, cm
import time
import numpy.ma as ma
from scipy.signal import savgol_filter
from matplotlib.ticker import AutoMinorLocator
import TM_Toolbox as TMT #personal Toolbox from Tobias Machnitzki
import matplotlib.patches as mpatches
import os
from PIL import Image
import glob
from matplotlib.collections import PatchCollection
#%%
#=========================
# # Parser:
#=========================

import argparse
parser = argparse.ArgumentParser()
 
parser.add_argument("-t", "--time", dest="yyyymmdd", type=str,
                    help="time of interest in the format yyyymmdd")
 
options = parser.parse_args()
 
yyyymmdd    = options.yyyymmdd
date_str    = yyyymmdd[2:8]
year_str    = yyyymmdd[0:4]
month_str   = yyyymmdd[4:6]
day_str     = yyyymmdd[6:8]

#%%
#=========================
# # Input NETCDF Files
#=========================

StartTime = 0
EndTime   = 24


nc_file = "/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/B_Reflectivity/Version_2/MMCR__MBR__Spectral_Moments__10s__155m-25km__"+date_str+".nc"
nc_meteo = "/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/I_Meteorology_2m/" + year_str + month_str + "/Meteorology__Deebles_Point__2m_10s__"+year_str + month_str + day_str +".nc" 
image_path = "/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_0/1_Allskyimager/data/cc" + year_str[2:4] + month_str + "/" + day_str + "/sm" + year_str[2:4] + month_str + day_str + ".tgz"

save_figures = "/home/mpim/m300517/MPI/working/Quicklook_live/save_figures/"
save_name = save_figures+'QL_Live'+year_str+month_str+day_str+'.png'

temp_folder = "temp_images_folder" #Name for a tempurary created folder for extracting the images

Filter = -55    #apply Filter of -55 dBz
MaxAlt = 20000  #set highest altitude to 20000 m

print('Load NetCDF Data')
print(nc_meteo)

#=========================
# # Read NETCDF variablespassing args to python script from python script
#=========================

#%%
#=========================
# # Get Data from Radar:
#=========================

fh = Dataset(nc_file, mode='r')

Z_dbz = fh.variables['Zf'][:].copy()	# in mm6 m-3
Time = fh.variables['time'][:].copy()
Range = fh.variables['range'][:].copy()
Status = fh.variables['status'][:].copy()

velocity = fh.variables['VEL'][:].copy()
print(fh.variables['time'])
fh.close()

# Umrechnung von Z in mm6 m-3 nach dBZ:
Z_dbz = np.where(ma.less_equal(Z_dbz,Filter),np.nan,Z_dbz)

# %%
# convert time from epoch to UTC for x-axis
secs = mdate.epoch2num(Time)

#%%
#================================
# # Get Data from meteorology 2m:
#================================

meteo = Dataset(nc_meteo)

temp = meteo.variables['T'][:].copy()
wind_speed = meteo.variables['VEL'][:].copy()
rain = meteo.variables['MXRI'][:].copy()
RH = meteo.variables['RH'][:].copy()
meteo_time = meteo.variables['time'][:].copy()

rain[rain < 0.05] = 0    #Filtering the error values

meteo_secs = mdate.epoch2num(meteo_time)

meteo.close()

# %%
#=========================
# # Get maxima and minima:
#=========================

print(Time[TMT.nanargmax(Z_dbz)[0]])
print(time.strftime('%H:%M:%S', time.gmtime((Time[TMT.nanargmax(Z_dbz)[0]])))[0:5])

MaxEcho = 'Highest Echo = ' +str(np.nanmax(Z_dbz))[0:5] + ' dBZ ('+ str(time.strftime('%H:%M:%S', time.gmtime((Time[TMT.nanargmax(Z_dbz)[0]])))[0:5]) + ' UTC)'

MinTemp = 'Min Temp = ' +str(np.nanmin(temp))[0:6] + ' °C ('+ str(time.strftime('%H:%M:%S', time.gmtime((meteo_time[TMT.nanargmin(temp)[0]])))[0:5]) + ' UTC)'
MaxTemp = 'Max Temp = ' +str(np.nanmax(temp))[0:6] + ' °C ('+ str(time.strftime('%H:%M:%S', time.gmtime((meteo_time[TMT.nanargmax(temp)[0]])))[0:5]) + ' UTC)'
DeltaTemp = '$\Delta$ Temp = ' + str(abs(abs(np.nanmax(temp)) - abs(np.nanmin(temp))))[0:6] + ' °C'
temp_str = MinTemp + '\n' + MaxTemp + '\n' + DeltaTemp
                                   
MaxRain = 'Max Precip = ' +str(np.nanmax(rain))[0:6] + ' mm/h ('+ str(time.strftime('%H:%M:%S', time.gmtime((meteo_time[TMT.nanargmax(rain)[0]])))[0:5]) + ' UTC)'


# %%
#=========================
# Setup Live-Data Strings:
#=========================
Humidity_test = 85.534
Precip_test = 24.123 

LiveTemp =  str(round(temp[-1],1)) + ' °C'
LivePrecip =  str(round(Precip_test,1)) + ' mm h$^{-1}$'
LiveHumid = str(round(Humidity_test,1)) + ' %'

live_time = mdate.num2date(meteo_secs[-1])
updatestr = live_time.strftime("%x %X")

wind_u= 30
wind_v = 20
#%%
#================================
# # apply Filter on reflectivity:
#================================

Z_dbz[Z_dbz < Filter] = np.nan
Z_dbz[Z_dbz > 40] = np.nan
plot_velocity = velocity.copy()
plot_velocity[plot_velocity > 3] = 3.5
plot_velocity[plot_velocity < -3] = -3.5
plot_velocity[np.isnan(Z_dbz)] = np.nan
     
#%%
#==========================================
# # Set starttime and endtime for one day:
#==========================================

StartTime=int((len(Time)/24*StartTime))
EndTime=int((len(Time)/24*EndTime))

#%%
#=======================================================
# Apply Savitzky-Golay filter on 2m-Data for smoothing:
#=======================================================

temp = savgol_filter(temp,61,1)
RH = savgol_filter(RH,61,1)
wind_speed = savgol_filter(wind_speed,61,1)

#%%
#===================================
# # Calculating Dewpointtemperature:
#===================================
tempd = TMT.Dewpoint(T=temp,RH=RH)



#%%
#===================================
# # Calculating LCL:
#===================================

z_0 = 0.02  #Hight of the surface measurement (station-hight in km)
LCL = (((temp-tempd)/8) + z_0 ) * 1000  # *1000 to convert from km to m


#%%
#===================================
# # Extracting Image:
#===================================
image_name_check = year_str[2:4] + month_str + day_str

if os.path.isdir(temp_folder):
    if len( glob.glob(temp_folder +"/*" + image_name_check + "*")) < 2:
        os.rmdir(temp_folder)
else:
    os.mkdir(temp_folder)
    
    tar_string = 'tar -xzvf ' + image_path  +' -C ' + temp_folder 
    os.system(tar_string)

wdir = os.getcwd()

image_dir = wdir + "/" + temp_folder + "/"
image_name = image_dir + sorted(os.listdir(image_dir))[720]

image = Image.open(image_name)

#%%
#===================================
# Making fancy boxes for Livedata
#===================================

def makeFancyBox(x,y,width,height):
    fancybox = mpatches.FancyBboxPatch(
    [x,y],width,height,
    boxstyle=mpatches.BoxStyle("Round", pad=1))
    patches.append(fancybox)

#%%
#====================
# Plot Figure
#====================

print('Plot Figure')

#Setting Fontsizes:
ax_title_size = 16
ylabel_size = 11
cb_title_size = 11
cb_size = 11
box_font = 11
legend_size = 11

barb_x = 73
barb_y = 16

#Building the frame and setting Title
fig = plt.figure(1,figsize=[16,9],facecolor="#CCCCCC")
gs1 = gridspec.GridSpec(56, 100)
plt.suptitle('Deebles Point, Barbados, '+day_str+'.'+month_str+'.'+year_str, fontsize=20,y=0.99)
plt.subplots_adjust(top=0.97,bottom=0.02,left =0.05,right=0.97)

patches=[]

ax1 = fig.add_subplot(gs1[5:20,:-40]) #Reflectivity
ax_cb = fig.add_subplot(gs1[5:20,50:62]) #Colorbar for Reflectivity
ax2 = fig.add_subplot(gs1[25:35,:-40], sharex = ax1) #Rain
ax3 = fig.add_subplot(gs1[40:50,:-40], sharex = ax1) #2m Temperature
#ax4 = ax3.twinx() #2m Relative Humidity
ax5 = ax3.twinx() #dewpoint temperature
ax6 = fig.add_subplot(gs1[25:50,71:94]) #Image
ax7 = fig.add_subplot(gs1[5:45,68:100]) #Live-Data frame layer
ax8 = fig.add_subplot(gs1[barb_y:barb_y+8,barb_x:barb_x+8]) #Live-Data-Windbarb
ax9 = fig.add_subplot(gs1[5:45,68:100],) #Live-Data Text layer
ax10 = fig.add_subplot(gs1[3:56,67:68]) #Black bar

#Adjusting appearance for all subplots:
axes = [ax1,ax2,ax3]
date_fmt = '%H:%M'
date_formatter = mdate.DateFormatter(date_fmt)

for ax in axes:
    plt.sca(ax)
    plt.xticks(rotation=70)

    ax.grid()

    ax.xaxis_date()
    ax.xaxis.set_major_formatter(date_formatter)
    ax.xaxis.set_minor_locator(AutoMinorLocator(6))

    ax.tick_params('both', length=3, width=1, which='minor',labelsize=ylabel_size)
    ax.tick_params('both', length=4, width=2, which='major',labelsize=ylabel_size)
    ax.get_yaxis().set_tick_params(which='both', direction='out')
    ax.get_xaxis().set_tick_params(which='both', direction='out')
    
ax3.set_xlabel('UTC Time',fontsize=ylabel_size)

#Setting up the x_axes for all plots:
ax3.set_xticks(np.arange(min(secs), max(secs),0.083333/2)) #Abstand xticks = 2h, 1 Minute = 0.000694
ax3.set_xlim(min(secs),max(secs))

#------------------------------------------------------------------------------------------------------------------------------------
#=========================
#Plotting Reflectivity:
#=========================   

print('Plotting Reflectivity...')
im1 = ax1.contourf(secs[StartTime:EndTime], Range[:], Z_dbz[StartTime:EndTime,:].transpose(), 100, vmin=Filter, vmax=40, cmap=cm.jet)
ax1.set_title('MBR2, Filter = '+str(Filter) + ' dBZ', fontsize=ax_title_size)


    #y-axes settings:
ax1.yaxis.set_ticks(np.arange(0,31000,2000))
ax1.set_ylim(0,MaxAlt)
ax1.set_ylabel('Altitude [km]',fontsize=ylabel_size)
ax1.set_yticklabels(np.arange(0,31,2))

    
#make x-axe-labels invisible:
for label in ax1.xaxis.get_ticklabels():
    label.set_visible(False)    

#Check if Radar-data is available everywhere:
LegendNoData = False
#plot vertical bars if the radar was off
for i in range(len(Time)):
    if Status[i] == 0:
        ax1.axvline(secs[i], linewidth = 1, color = 'grey')
        LegendNoData = True

#Plot LCL:
lcl_im, = ax1.plot(meteo_secs,LCL,color='black', ls = '--', lw=1.5, label = 'Lifting Condensation Level')
lcl_legend = ax1.legend(loc=1, fontsize=legend_size)

if LegendNoData:
    grey_patch = mpatches.Patch(color='grey', label='RADAR off')
    grey_legend = ax1.legend(handles=[lcl_im,grey_patch], fontsize=legend_size, loc=1)

    #collorbar settings:
ax_cb.set_xticks([])
ax_cb.set_yticks([])
ax_cb.set_visible(False)
cb1 = plt.colorbar(im1, ax = ax_cb, ticks=[x for x in range(-50,41,10)], shrink=1)
cb1.set_clim(Filter, 40)
cb1.set_label('Radar reflectivity factor Zf [dBZ]',fontsize=cb_title_size)
cb1.ax.tick_params(labelsize=cb_size)

    #Creating box with Max and Min values:
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(0.67, 1.1, MaxEcho, transform=ax1.transAxes, fontsize=box_font, verticalalignment='top', bbox=props)


#------------------------------------------------------------------------------------------------------------------------------------
#=========================
#Plotting Rain:
#=========================   
print('Plotting Rain...')
Legend_handle = False
for i in range(len(meteo_secs)):
    if rain[i] > 0.1:
        ax2.axvline(meteo_secs[i], linewidth = 1, color = 'g',alpha=0.025, label='Timeframe where precipitation was measured')
        Legend_handle = True  

im2, = ax2.plot(meteo_secs,rain,color='b',lw=1, label='Precipitation')
ax2.set_title('Precipitation', fontsize=ax_title_size)

#make x-axe-labels invisible:
for label in ax2.xaxis.get_ticklabels():
    label.set_visible(False) 

    #y-axes settings:
ax2.set_ylabel('Precipitation [mm/h]', fontsize=ylabel_size)
ax2.set_ylim(0,50)
ax2.set_yticks(np.arange(0,51,5))
ax2.legend(loc='upper right', fontsize = legend_size)

    #make shaded area where Rain was measured:        
if Legend_handle == True:
    green_patch = mpatches.Patch(color='g',alpha = 0.2, label='Timeframe where precipitation was measured')
    green_legend = ax2.legend(handles=[im2,green_patch], fontsize=legend_size, loc=1)     

if np.nanmax(rain) > 50:
    rain_max = round(np.nanmax(rain)+5 , 0)
    ax2.set_ylim(0,rain_max)
    ax2.set_yticks(np.arange(0,rain_max,5))
    for label in ax2.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False) #making everey 2nd label invisibel 


    #Creating box with Max and Min values:
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(0.67, 1.15, MaxRain, transform=ax2.transAxes, fontsize=box_font, verticalalignment='top', bbox=props)


#------------------------------------------------------------------------------------------------------------------------------------
#======================================
#Plotting Temperature, Dewpoint and RH:
#======================================
'''
Temperature, Dewpoint and RH are plotted all in seperate axes, but in one subplot.
This means, that they all share the same x-axes.
Furthermore The Temperature and Dewpoint share the same y-axes (left) and
the Relative Humidity is on the right y-axes.
'''
print('Plotting Temperature...')

#----------------------------
#Temperature
#----------------------------
im3 = ax3.plot(meteo_secs,temp, color='red', lw=1, label='Temperature')
ax3.set_title('2 m Temperature', fontsize=ax_title_size)

    #y-axes settings:
ax3.yaxis.set_ticks(np.arange(-20,70,1))    
ax3.set_ylim(min(tempd)-1,max(temp)+1)
ax3.set_ylabel('Temperature [°C]', fontsize=ylabel_size)

    #creating a grid for the Temperature and Dewpointtemperature:
#ax3.grid(color = 'grey',lw=1)
    #Creating box with Max and Min values:
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(0.67, 1.15, MaxTemp, transform=ax3.transAxes, fontsize=box_font, verticalalignment='top', bbox=props)


#----------------------------
#Dewpointtemperature:
#----------------------------
im5=ax5.plot(meteo_secs,tempd,color='b',lw=1, label='Dewpoint temperature')
ax5.yaxis.tick_left()
    #y-axes settings, so that Temperature and Dewpointtemperature share the same y-axes:
ax5.set_ylim(min(tempd)-1,max(temp)+1)
for label in ax5.yaxis.get_ticklabels():
    label.set_visible(False)    #making the y-labels invisbel as they are the same as the temperature-labels


##----------------------------
##Relative Humidity
##----------------------------
#im4 = ax4.plot(meteo_secs,RH,color='g',lw=1, label='Relative Humidity [%]')
#
#    #creating shaded area between RH-line and x-axes
#ax4.fill_between(meteo_secs,0,RH,facecolor='g',alpha=0.2)
#")
#    #y-axes settings. The y-axes for RH is on the right side!:
#ax4.set_ylabel('Relative Humidity [%]', fontsize= ylabel_size)
#ax4.yaxis.set_ticks(np.arange(0,101,10))
#ax4.set_ylim(0, 110)
#ax4.tick_params('y', length=4, width=2, which='major',labelsize=ylabel_size)
#ax4.yaxis.set_label_coords(1.04,0.5)
#
#    #creating a dashed-grid for the Relative Humidity:
#ax4.grid(color = 'grey',lw=1,ls='--', axis='y')
#
##------------------------------------------------
##Creating the Legend for the Temp. Dewp. and RH:
##------------------------------------------------
#lns = im3+im5+im4   #putting together all 3 axes in one list
#labs = [l.get_label() for l in lns] #Getting the label from each plot and put it in a list
#ax5.legend(lns, labs, fontsize=legend_size, loc='lower left', borderaxespad=0.) #creating one legend for the 3 axes (ax3, ax4, ax8)

lns = im3+im5   #putting together all 3 axes in one list
labs = [l.get_label() for l in lns] #Getting the label from each plot and put it in a list
ax5.legend(lns, labs, fontsize=legend_size, loc=1) #creating one legend for the 3 axes (ax3, ax4, ax8)


#======================================
#Plotting The Image:
#======================================
print('Plotting Image...')
ax6.imshow(image)
ax6.set_xticks([])
ax6.set_yticks([])



#======================================
# Live-Data:
#======================================
print('Plotting Live-Data...')

live_number_font = box_font +5

ax7.axis('off') #Creating just a field for text to be drawn into
ax9.axis('off')

ax7.set_zorder(7)
ax8.set_zorder(9)
ax9.set_zorder(9)

#ax7.axis('equal')
ax7.set_xlim(0,100)
ax7.set_ylim(0,100)

ax9.set_xlim(0,100)
ax9.set_ylim(0,100)



makeFancyBox(10,93,71,5)  #Title
makeFancyBox(10,74,33,15)  # Temperature
makeFancyBox(48,74,33,15)  #Precipitation
makeFancyBox(10,55,33,15)  # Wind
makeFancyBox(48,55,33,15)  #Relative Humidity


collection = PatchCollection(patches,color='white',edgecolor='black')
ax7.add_collection(collection)


ax9.text(0, 0.97, '                    Live Data            ', transform=ax9.transAxes, fontsize=box_font+7, verticalalignment='top')
ax9.text(0, 0.98, ('last update ' + updatestr), transform=ax9.transAxes, fontsize=box_font+7, verticalalignment='top')


ax9.text(0.1, 0.89,("Temperature"), transform=ax9.transAxes, fontsize=box_font, verticalalignment='top',color='black')
ax9.text(0.17, 0.83, (LiveTemp), transform=ax9.transAxes, fontsize=live_number_font, verticalalignment='top')

ax9.text(0.48, 0.89,("Precipitation"), transform=ax9.transAxes, fontsize=box_font, verticalalignment='top')
if Precip_test < 100:   #if precipitation > 99.9: make the fontsize smaller to fit in the box
    precip_font = 0
else:
    precip_font = -1    
ax9.text(0.5, 0.83, (LivePrecip), transform=ax9.transAxes, fontsize=(live_number_font + precip_font), verticalalignment='top')


ax9.text(0.48, 0.7,("Relative Humidity"), transform=ax9.transAxes, fontsize=box_font, verticalalignment='top')
ax9.text(0.56, 0.63, (LiveHumid), transform=ax9.transAxes, fontsize=live_number_font, verticalalignment='top')

    

ax9.text(0.1, 0.7,("Wind"), transform=ax9.transAxes, fontsize=box_font, verticalalignment='top')

ax8.set_visible(True)
ax8.axis('off')
ax8.barbs(1,1,wind_u,wind_v,pivot='middle', length=10) # wind_u and wind_v must be in knots
ax8.set_xlim(0,2)
ax8.set_ylim(0,2)
ax8.set_facecolor('white')

#======================================
#Black bar for seperation:
#======================================
#ax10.axis('off')
ax10.get_xaxis().set_visible(False)
ax10.get_yaxis().set_visible(False)
ax10.set_zorder(10)
ax10.set_facecolor('black')
ax10.plot(1)

#%%
#======================================
#Saving the Figure:
#======================================
print('saving figure...')
plt.savefig(save_name, facecolor=fig.get_facecolor(), edgecolor='none')
print('saved')

