#!/usr/bin/python
# -*- coding: latin-1 -*-
import numpy as np
import scipy as sci
import math as m
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate

#Define function to return solar flux for a given day of the year
def solarflux(n):
    solar_const=1367.0
    return solar_const*(1+0.0034*(np.cos(2*m.pi*(n-1)/365.25)))

#Define function to return k and c values for a given month of the year
def ck(m):
    k=[0.142,0.144,0.156,0.18,0.196,0.205,0.207,0.201,0.177,0.160,0.149,0.142]
    c=[0.058,0.06,0.071,0.097,0.121,0.134,0.136,0.122,0.092,0.073,0.063,0.057]
    return np.exp(-k[m-1]),c[m-1]

#Read data from CSV file and store into array
solar_data_night=np.genfromtxt("NSRDB_StationData_20020101_20021231_723870.csv",delimiter=',',skip_header=1)

#Calculate total number of rows in this array that correspond to daytime
count_day=0
for i in range(len(solar_data_night)):
    if solar_data_night[i,2]!=99:
        count_day=count_day+1

#Create zero array to store data for daytimes
solar_data_day=np.zeros(shape=(count_day,13),dtype=float)

#Store daytime zenith and azimuthal angles into this array
count=0
for i in range(len(solar_data_night)):
    if solar_data_night[i,2]!=99:
        solar_data_day[count,0]=solar_data_night[i,0]
        solar_data_day[count,1]=solar_data_night[i,1]
        solar_data_day[count,2]=solar_data_night[i,2]
        solar_data_day[count,3]=solar_data_night[i,3]
        count=count+1
    else:
        continue

#Initialise day counter
day=1
    
#Calculate and store solar radiation parameters into the array
for i in range(len(solar_data_day)):
    
    #Calculate month
    month=int(solar_data_day[i,0]/100)
    
    #Calculate zenith angle in radians    
    zen_ang=solar_data_day[i,2]*(m.pi/180.0)
    
    #Calculate azimuthal angle in radians
    azi_ang=(180-solar_data_day[i,3])*(m.pi/180.0)
    
    #Calculate surface azimuthal angle in radians    
    surfazi=0*(m.pi/180.0)
    
    #Calculate surface zenith angle in radians
    surftheta=40*(m.pi/180.0)
    
    #Chose reflectivity, based on albedo values in http://en.wikipedia.org/wiki/Albedo
    #We have chosen rho=0.4 for desert sand
    rho=0.4
    
    #Calculate day
    if(i==0):
        x=101
    elif(i>0):
        x=solar_data_day[i-1,0]
    if (solar_data_day[i,0]-x>0): 
        day=day+1    
    
    #Calculate and store Gbn
    solar_data_day[i,4]=solarflux(day)*(ck(month)[0]**(1.0/np.cos(zen_ang)))
    #Calculate and store Gbh
    solar_data_day[i,5]=solar_data_day[i,4]*np.cos(zen_ang)
    #Calculate and store Gdh
    solar_data_day[i,6]=solar_data_day[i,4]*ck(month)[1]
    #Calculate and store Gh
    solar_data_day[i,7]=solar_data_day[i,5]+solar_data_day[i,6]
    #Calculate and store angle of incidence
    solar_data_day[i,8]=np.arccos(np.cos(zen_ang)*np.cos(surftheta)+np.sin(zen_ang)*np.sin(surftheta)*np.cos(azi_ang-surfazi))
    #Calculate and store Gbt    
    solar_data_day[i,9]=solar_data_day[i,4]*np.cos(solar_data_day[i,8])
    #Calculate and store Gdt
    solar_data_day[i,10]=solar_data_day[i,4]*ck(month)[1]*0.5*(1+np.cos(surftheta))
    #Calculate and store Grt
    solar_data_day[i,11]=solar_data_day[i,7]*rho*0.5*(1-np.cos(surftheta))
    #Calculate and store Gt
    solar_data_day[i,12]=solar_data_day[i,9]+solar_data_day[i,10]+solar_data_day[i,11]

#Save hourly radiation data to file
np.savetxt("Hourly_Data.csv",solar_data_day,delimiter=',',header="Date,Time,Zenith Angle,Azi Angle,Gbn,Gbh,Gdh,Gh,i,Gbt,Gdt,Grt,Gt")

#Create new array to store daily incident radiation
daily_data=np.zeros(shape=(365,2),dtype=float)

#Calculate daily radiation incident upon collector
count=0
for i in range(len(solar_data_day)):
    if(i==0):
        x=101
    elif(i>0):
        x=solar_data_day[i-1,0]
    
    if(solar_data_day[i,0]-x==0):
        daily_data[count,0]=count+1
        daily_data[count,1]=daily_data[count,1]+solar_data_day[i,12]
    else:
        count=count+1
        daily_data[count,0]=solar_data_day[i,0]
        daily_data[count,1]=daily_data[count,1]+solar_data_day[i,12]

#Save daily radiation data to file
np.savetxt("Daily_Data.csv",daily_data,delimiter=',',header="Day,Incident Radiation",fmt='%3f')
        
#Create interpolating function, used later in plot formatting
f=interpolate.interp1d(daily_data[:,0],daily_data[:,1],kind='cubic')

#Plot
fig=plt.figure()
fig.suptitle("Daily Variation of Incident Solar Radiation")
ax1=fig.add_subplot(111)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Day of the Year")
ax1.set_ylim([6000,10000])
ax1.set_xlim([1,365])
ax1.set_title("On Collector Inclined at 40deg, facing South")
ax1.plot(daily_data[:,0],daily_data[:,1],color='r',linewidth=2)
ax1.text(79, 9500, "Incident Rad @ Spr. Equinox = %.2f" % f(79),rotation='vertical')
ax1.plot([79,79],[6000,10000],'k-')
ax1.text(172, 9500, "Incident Rad @ Sum. Solstice = %.2f" % f(172),rotation='vertical')
ax1.plot([172,172],[6000,10000],'k-')
ax1.text(265, 9500, "Incident Rad @ Fall Equinox = %.2f" % f(265),rotation='vertical')
ax1.plot([265,265],[6000,10000],'k-')
ax1.text(355, 9500, "Incident Rad @ Win. Equinox = %.2f" % f(355),rotation='vertical')
ax1.plot([355,355],[6000,10000],'k-')
fig.savefig("Daily Variation of Incident Solar Radiation.png")


    
    
    
    