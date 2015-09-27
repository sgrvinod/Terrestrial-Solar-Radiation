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

#Calculate number of columns
count_col=len(solar_data_night[0])

#Create zero array to store data for daytimes
solar_data_day=np.zeros(shape=(count_day,count_col),dtype=float)

#Specify 5 days of interest
a=104
b=323
c=622
d=918
e=1222

#Initialise count for calculating number of daylight hours in each of the 5 days
counthours=[0,0,0,0,0]

#Store daytime data into this new array
count=0
for i in range(len(solar_data_night)):
    if solar_data_night[i,2]!=99:
        solar_data_day[count]=solar_data_night[i]
        
        #Also calculate number of daylight hours in each of the 5 days
        if(solar_data_day[count,0]==a):
            counthours[0]+=1 
        if(solar_data_day[count,0]==b):
            counthours[1]+=1 
        if(solar_data_day[count,0]==c):
            counthours[2]+=1 
        if(solar_data_day[count,0]==d):
            counthours[3]+=1 
        if(solar_data_day[count,0]==e):
            counthours[4]+=1 
        count=count+1
    
    else:
        continue

#Create new array to store HOURLY measured, modeled and c-k modeled data
all_data=np.zeros(shape=(4820,13),dtype=float)

#Create new array to store DAILY measured, modeled and c-k modeled data for the 5 specified days
day_points=np.zeros(shape=(5,10),dtype=float)

#Create new array to store error in NSRDB and C-k models
error=np.zeros(shape=(5,7),dtype=float)

#Create new arrays to store HOURLY values for the 5 days of interest
day1=np.zeros(shape=(counthours[0],13),dtype=float)
day2=np.zeros(shape=(counthours[1],13),dtype=float)
day3=np.zeros(shape=(counthours[2],13),dtype=float)
day4=np.zeros(shape=(counthours[3],13),dtype=float)
day5=np.zeros(shape=(counthours[4],13),dtype=float)

#Initialise hour count for the 5 days, used later for assignment
counth=0

#Initialise day counter
day=1

#Calculate HOURLY measured, modeled and c-k modeled data
for i in range(len(all_data)):
    
    #Calculate month
    month=int(solar_data_day[i,0]/100)
    
    #Calculate zenith angle in radians    
    zen_ang=solar_data_day[i,2]*(m.pi/180.0)
    
    #Calculate day
    if(i==0):
        x=101
    elif(i>0):
        x=solar_data_day[i-1,0]
    if (solar_data_day[i,0]-x>0): 
        day=day+1
    
    #Store date, HOURLY time, zenith angle, azimuthal angle
    all_data[i,0]=solar_data_day[i,0]
    all_data[i,1]=solar_data_day[i,1]
    all_data[i,2]=solar_data_day[i,2]
    all_data[i,3]=solar_data_day[i,3]
    
    #Store HOURLY measured normal, diffuse and total radiation
    all_data[i,4]=solar_data_day[i,17]
    all_data[i,5]=solar_data_day[i,19]
    all_data[i,6]=solar_data_day[i,15]
    
    #Store HOURLY modeled normal, diffuse and total radiation
    all_data[i,7]=solar_data_day[i,9]
    all_data[i,8]=solar_data_day[i,12]
    all_data[i,9]=solar_data_day[i,6]
    
    #Calculate HOURLY c-k modeled normal, diffuse and total radiation
    all_data[i,10]=solarflux(day)*(ck(month)[0]**(1.0/np.cos(zen_ang)))
    all_data[i,11]=all_data[i,10]*ck(month)[1]
    all_data[i,12]=all_data[i,11]+all_data[i,10]*np.cos(zen_ang)
    
    #Calculate DAILY measured, modeled, c-k modeled data for 1st day of interest
    if(solar_data_day[i,0]==a):
        day_points[0,0]=a
        
        #Also store HOURLY data for 1st day of interest
        day1[counth]=all_data[i]
        counth+=1
        
        #Measured
        day_points[0,1]=day_points[0,1]+all_data[i,4]
        day_points[0,2]=day_points[0,2]+all_data[i,5]
        day_points[0,3]=day_points[0,3]+all_data[i,6]
        #Modeled        
        day_points[0,4]=day_points[0,4]+all_data[i,7]
        day_points[0,5]=day_points[0,5]+all_data[i,8]
        day_points[0,6]=day_points[0,6]+all_data[i,9]
        #c-k modeled        
        day_points[0,7]=day_points[0,7]+all_data[i,10]
        day_points[0,8]=day_points[0,8]+all_data[i,11]
        day_points[0,9]=day_points[0,9]+all_data[i,12]
           
    #Calculate DAILY measured, modeled, c-k modeled data for 2nd day of interest
    if(solar_data_day[i,0]==b):
        day_points[1,0]=b
        
        #Also store HOURLY data for 2nd day of interest
        day2[counth-counthours[0]]=all_data[i]
        counth+=1
        
        #Measured
        day_points[1,1]=day_points[1,1]+all_data[i,4]
        day_points[1,2]=day_points[1,2]+all_data[i,5]
        day_points[1,3]=day_points[1,3]+all_data[i,6]
        #Modeled
        day_points[1,4]=day_points[1,4]+all_data[i,7]
        day_points[1,5]=day_points[1,5]+all_data[i,8]
        day_points[1,6]=day_points[1,6]+all_data[i,9]
        #c-k modeled
        day_points[1,7]=day_points[1,7]+all_data[i,10]
        day_points[1,8]=day_points[1,8]+all_data[i,11]
        day_points[1,9]=day_points[1,9]+all_data[i,12]
        
    #Calculate DAILY measured, modeled, c-k modeled data for 3rd day of interest
    if(solar_data_day[i,0]==c):
        day_points[2,0]=c
        
        #Also store HOURLY data for 3rd day of interest
        day3[counth-counthours[0]-counthours[1]]=all_data[i]
        counth+=1
        
        #Measured
        day_points[2,1]=day_points[2,1]+all_data[i,4]
        day_points[2,2]=day_points[2,2]+all_data[i,5]
        day_points[2,3]=day_points[2,3]+all_data[i,6]
        #Modeled
        day_points[2,4]=day_points[2,4]+all_data[i,7]
        day_points[2,5]=day_points[2,5]+all_data[i,8]
        day_points[2,6]=day_points[2,6]+all_data[i,9]
        #c-k modeled
        day_points[2,7]=day_points[2,7]+all_data[i,10]
        day_points[2,8]=day_points[2,8]+all_data[i,11]
        day_points[2,9]=day_points[2,9]+all_data[i,12]
        
    #Calculate DAILY measured, modeled, c-k modeled data for 4th day of interest
    if(solar_data_day[i,0]==d):
        day_points[3,0]=d
        
        #Also store HOURLY data for 4th day of interest
        day4[counth-counthours[0]-counthours[1]-counthours[2]]=all_data[i]
        counth+=1
        
        #Measured
        day_points[3,1]=day_points[3,1]+all_data[i,4]
        day_points[3,2]=day_points[3,2]+all_data[i,5]
        day_points[3,3]=day_points[3,3]+all_data[i,6]
        #Modeled
        day_points[3,4]=day_points[3,4]+all_data[i,7]
        day_points[3,5]=day_points[3,5]+all_data[i,8]
        day_points[3,6]=day_points[3,6]+all_data[i,9]
        #c-k modeled        
        day_points[3,7]=day_points[3,7]+all_data[i,10]
        day_points[3,8]=day_points[3,8]+all_data[i,11]
        day_points[3,9]=day_points[3,9]+all_data[i,12]
        
    #Calculate DAILY measured, modeled, c-k modeled data for 5th day of interest
    if(solar_data_day[i,0]==e):
        day_points[4,0]=e
        
        #Also store HOURLY data for 5th day of interest
        day5[counth-counthours[0]-counthours[1]-counthours[2]-counthours[3]]=all_data[i]
        counth+=1
        
        #Measured
        day_points[4,1]=day_points[4,1]+all_data[i,4]
        day_points[4,2]=day_points[4,2]+all_data[i,5]
        day_points[4,3]=day_points[4,3]+all_data[i,6]
        #Modeled        
        day_points[4,4]=day_points[4,4]+all_data[i,7]
        day_points[4,5]=day_points[4,5]+all_data[i,8]
        day_points[4,6]=day_points[4,6]+all_data[i,9]
        #c-k modeled        
        day_points[4,7]=day_points[4,7]+all_data[i,10]
        day_points[4,8]=day_points[4,8]+all_data[i,11]
        day_points[4,9]=day_points[4,9]+all_data[i,12]
        
#Calculate error
for i in range(5):
    error[i,0]=day_points[i,0]
    error[i,1]=100*(day_points[i,4]-day_points[i,1])/day_points[i,1]
    error[i,2]=100*(day_points[i,5]-day_points[i,2])/day_points[i,2]
    error[i,3]=100*(day_points[i,6]-day_points[i,3])/day_points[i,3]
    error[i,4]=100*(day_points[i,7]-day_points[i,1])/day_points[i,1]
    error[i,5]=100*(day_points[i,8]-day_points[i,2])/day_points[i,2]
    error[i,6]=100*(day_points[i,9]-day_points[i,3])/day_points[i,3]
                
print "Errors (in %) are also written to a CSV file."
print error
             
#Save HOURLY measured, modeled and c-k modeled data
np.savetxt("Hourly_Data.csv",all_data,delimiter=',',header="Date,Time,Zenith Angle,Azimuthal Angle,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")

#Save HOURLY measured, modeled and c-k modeled data for the 5 days of interest
np.savetxt("Hourly_Day1.csv",day1,delimiter=',',header="Date,Time,Zenith Angle,Azimuthal Angle,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")
np.savetxt("Hourly_Day2.csv",day2,delimiter=',',header="Date,Time,Zenith Angle,Azimuthal Angle,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")
np.savetxt("Hourly_Day3.csv",day3,delimiter=',',header="Date,Time,Zenith Angle,Azimuthal Angle,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")
np.savetxt("Hourly_Day4.csv",day4,delimiter=',',header="Date,Time,Zenith Angle,Azimuthal Angle,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")
np.savetxt("Hourly_Day5.csv",day5,delimiter=',',header="Date,Time,Zenith Angle,Azimuthal Angle,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")

#Save DAILY measured, modeled and c-k modeled data for specified days
np.savetxt("Daily_Data.csv",day_points,delimiter=',',header="Date,Meas Normal,Meas Diff,Meas Glo,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo")

#Save model error
np.savetxt("Models'_Error_in%.csv",error,delimiter=',',header="Date,Mod Normal,Mod Diff,Mod Glo,Ck Normal,Ck Diff,Ck Glo",fmt="%4f")

#Plot for day 1
fig=plt.figure()
fig.suptitle("A COMPARISION for %d (MMDD) (Red - Measured | Blue - NSRDB Model | Green - ck Model)" % a)
fig.text(0.5,0.75,"Discrepancies in the red and blue curves\nw.r.t green are caused by varying cloud cover",ha='center',style='italic',fontsize='small')

ax1=fig.add_subplot(311)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Hour")
ax1.set_ylim([0,1100])
ax1.set_title("Normal")
ax1.plot(day1[:,1],day1[:,4],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax5=plt.twinx(ax1)
ax5.set_ylim([0,1100])
ax5.plot(day1[:,1],day1[:,7],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax6=plt.twinx(ax1)
ax6.set_ylim([0,1100])
ax6.plot(day1[:,1],day1[:,10],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax2=fig.add_subplot(312)
ax2.set_ylabel("Wh/m2")
ax2.set_xlabel("Hour")
ax2.set_ylim([0,350])
ax2.set_title("Diffused")
ax2.plot(day1[:,1],day1[:,5],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax7=plt.twinx(ax2)
ax7.set_ylim([0,350])
ax7.plot(day1[:,1],day1[:,8],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax8=plt.twinx(ax2)
ax8.set_ylim([0,350])
ax8.plot(day1[:,1],day1[:,11],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax3=fig.add_subplot(313)
ax3.set_ylabel("Wh/m2")
ax3.set_xlabel("Hour")
ax3.set_ylim([0,650])
ax3.set_title("Global")
ax3.plot(day1[:,1],day1[:,6],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax9=plt.twinx(ax3)
ax9.set_ylim([0,650])
ax9.plot(day1[:,1],day1[:,9],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax10=plt.twinx(ax3)
ax10.set_ylim([0,650])
ax10.plot(day1[:,1],day1[:,12],color='g',label="Emissivity vs. Wavelength",linewidth=2)

fig.subplots_adjust(wspace=0.6, hspace=0.7)
fig.savefig("Solar_Radiation_Model_Comparision_for_%d.png" % a)

#Plot for day 2
fig=plt.figure()
fig.suptitle("A COMPARISION for %d (MMDD) (Red - Measured | Blue - NSRDB Model | Green - ck Model)" % b)
fig.text(0.5,0.75,"Discrepancies in the red and blue curves\nw.r.t green are caused by varying cloud cover",ha='center',style='italic',fontsize='small')

ax1=fig.add_subplot(311)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Hour")
ax1.set_ylim([0,1200])
ax1.set_title("Normal")
ax1.plot(day2[:,1],day2[:,4],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax5=plt.twinx(ax1)
ax5.set_ylim([0,1200])
ax5.plot(day2[:,1],day2[:,7],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax6=plt.twinx(ax1)
ax6.set_ylim([0,1200])
ax6.plot(day2[:,1],day2[:,10],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax2=fig.add_subplot(312)
ax2.set_ylabel("Wh/m2")
ax2.set_xlabel("Hour")
ax2.set_ylim([0,350])
ax2.set_title("Diffused")
ax2.plot(day2[:,1],day2[:,5],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax7=plt.twinx(ax2)
ax7.set_ylim([0,350])
ax7.plot(day2[:,1],day2[:,8],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax8=plt.twinx(ax2)
ax8.set_ylim([0,350])
ax8.plot(day2[:,1],day2[:,11],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax3=fig.add_subplot(313)
ax3.set_ylabel("Wh/m2")
ax3.set_xlabel("Hour")
ax3.set_ylim([0,1100])
ax3.set_title("Global")
ax3.plot(day2[:,1],day2[:,6],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax9=plt.twinx(ax3)
ax9.set_ylim([0,1100])
ax9.plot(day2[:,1],day2[:,9],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax10=plt.twinx(ax3)
ax10.set_ylim([0,1100])
ax10.plot(day2[:,1],day2[:,12],color='g',label="Emissivity vs. Wavelength",linewidth=2)

fig.subplots_adjust(wspace=0.6, hspace=0.7)
fig.savefig("Solar_Radiation_Model_Comparision_for_%d.png" % b)

#Plot for day 3
fig=plt.figure()
fig.suptitle("A COMPARISION for %d (MMDD) (Red - Measured | Blue - NSRDB Model | Green - ck Model)" % c)
fig.text(0.5,0.75,"Discrepancies in the red and blue curves\nw.r.t green are caused by varying cloud cover",ha='center',style='italic',fontsize='small')

ax1=fig.add_subplot(311)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Hour")
ax1.set_ylim([0,1200])
ax1.set_title("Normal")
ax1.plot(day3[:,1],day3[:,4],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax5=plt.twinx(ax1)
ax5.set_ylim([0,1200])
ax5.plot(day3[:,1],day3[:,7],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax6=plt.twinx(ax1)
ax6.set_ylim([0,1200])
ax6.plot(day3[:,1],day3[:,10],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax2=fig.add_subplot(312)
ax2.set_ylabel("Wh/m2")
ax2.set_xlabel("Hour")
ax2.set_ylim([0,350])
ax2.set_title("Diffused")
ax2.plot(day3[:,1],day3[:,5],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax7=plt.twinx(ax2)
ax7.set_ylim([0,350])
ax7.plot(day3[:,1],day3[:,8],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax8=plt.twinx(ax2)
ax8.set_ylim([0,350])
ax8.plot(day3[:,1],day3[:,11],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax3=fig.add_subplot(313)
ax3.set_ylabel("Wh/m2")
ax3.set_xlabel("Hour")
ax3.set_ylim([0,1300])
ax3.set_title("Global")
ax3.plot(day3[:,1],day3[:,6],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax9=plt.twinx(ax3)
ax9.set_ylim([0,1300])
ax9.plot(day3[:,1],day3[:,9],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax10=plt.twinx(ax3)
ax10.set_ylim([0,1300])
ax10.plot(day3[:,1],day3[:,12],color='g',label="Emissivity vs. Wavelength",linewidth=2)

fig.subplots_adjust(wspace=0.6, hspace=0.7)
fig.savefig("Solar_Radiation_Model_Comparision_for_%d.png" % c)

#Plot for day 4
fig=plt.figure()
fig.suptitle("A COMPARISION for %d (MMDD) (Red - Measured | Blue - NSRDB Model | Green - ck Model)" % d)
fig.text(0.5,0.75,"Discrepancies in the red and blue curves\nw.r.t green are caused by varying cloud cover",ha='center',style='italic',fontsize='small')

ax1=fig.add_subplot(311)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Hour")
ax1.set_ylim([0,1200])
ax1.set_title("Normal")
ax1.plot(day4[:,1],day4[:,4],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax5=plt.twinx(ax1)
ax5.set_ylim([0,1200])
ax5.plot(day4[:,1],day4[:,7],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax6=plt.twinx(ax1)
ax6.set_ylim([0,1200])
ax6.plot(day4[:,1],day4[:,10],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax2=fig.add_subplot(312)
ax2.set_ylabel("Wh/m2")
ax2.set_xlabel("Hour")
ax2.set_ylim([0,350])
ax2.set_title("Diffused")
ax2.plot(day4[:,1],day4[:,5],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax7=plt.twinx(ax2)
ax7.set_ylim([0,350])
ax7.plot(day4[:,1],day4[:,8],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax8=plt.twinx(ax2)
ax8.set_ylim([0,350])
ax8.plot(day4[:,1],day4[:,11],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax3=fig.add_subplot(313)
ax3.set_ylabel("Wh/m2")
ax3.set_xlabel("Hour")
ax3.set_ylim([0,1100])
ax3.set_title("Global")
ax3.plot(day4[:,1],day4[:,6],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax9=plt.twinx(ax3)
ax9.set_ylim([0,1100])
ax9.plot(day4[:,1],day4[:,9],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax10=plt.twinx(ax3)
ax10.set_ylim([0,1100])
ax10.plot(day4[:,1],day4[:,12],color='g',label="Emissivity vs. Wavelength",linewidth=2)

fig.subplots_adjust(wspace=0.6, hspace=0.7)
fig.savefig("Solar_Radiation_Model_Comparision_for_%d.png" % d)

#Plot for day 5
fig=plt.figure()
fig.suptitle("A COMPARISION for %d (MMDD) (Red - Measured | Blue - NSRDB Model | Green - ck Model)" % e)
fig.text(0.5,0.75,"Discrepancies in the red and blue curves\nw.r.t green are caused by varying cloud cover",ha='center',style='italic',fontsize='small')

ax1=fig.add_subplot(311)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Hour")
ax1.set_ylim([0,1200])
ax1.set_title("Normal")
ax1.plot(day5[:,1],day5[:,4],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax5=plt.twinx(ax1)
ax5.set_ylim([0,1200])
ax5.plot(day5[:,1],day5[:,7],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax6=plt.twinx(ax1)
ax6.set_ylim([0,1200])
ax6.plot(day5[:,1],day5[:,10],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax2=fig.add_subplot(312)
ax2.set_ylabel("Wh/m2")
ax2.set_xlabel("Hour")
ax2.set_ylim([0,350])
ax2.set_title("Diffused")
ax2.plot(day5[:,1],day5[:,5],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax7=plt.twinx(ax2)
ax7.set_ylim([0,350])
ax7.plot(day5[:,1],day5[:,8],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax8=plt.twinx(ax2)
ax8.set_ylim([0,350])
ax8.plot(day5[:,1],day5[:,11],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax3=fig.add_subplot(313)
ax3.set_ylabel("Wh/m2")
ax3.set_xlabel("Hour")
ax3.set_ylim([0,650])
ax3.set_title("Global")
ax3.plot(day5[:,1],day5[:,6],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax9=plt.twinx(ax3)
ax9.set_ylim([0,650])
ax9.plot(day5[:,1],day5[:,9],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax10=plt.twinx(ax3)
ax10.set_ylim([0,650])
ax10.plot(day5[:,1],day5[:,12],color='g',label="Emissivity vs. Wavelength",linewidth=2)

fig.subplots_adjust(wspace=0.6, hspace=0.7)
fig.savefig("Solar_Radiation_Model_Comparision_for_%d.png" % e)

#Plot for the year 2002
fig=plt.figure()
fig.suptitle("A PAN-YEAR COMPARISION (Red - Measured | Blue - NSRDB Model | Green - ck Model)")
fig.text(0.5,0.73,"We have chosen 1/04 (CloudCover~0.01), 3/23 (CC~0.35),\n6/22 (CC~0), 12/23 (CC~0.6) instead of 1/01 (CC~0.52), 3/21 (CC~0.88),\n6/21 (CC~0.59), 12/21 (CC~0.75)",ha='center',style='italic',fontsize='small')

ax1=fig.add_subplot(311)
ax1.set_ylabel("Wh/m2")
ax1.set_xlabel("Day of the Year")
ax1.set_ylim([0,15000])
ax1.set_title("Normal")
x=[1,2,3,4,5]
my_xticks=["01/03","03/23","06/22","09/18","12/22"]
ax1.set_xticks(x)
ax1.set_xticklabels(my_xticks)
ax1.plot(x,day_points[:,1],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax5=plt.twinx(ax1)
ax5.set_ylim([0,15000])
ax5.plot(x,day_points[:,4],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax6=plt.twinx(ax1)
ax6.set_ylim([0,15000])
ax6.plot(x,day_points[:,7],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax2=fig.add_subplot(312)
ax2.set_ylabel("Wh/m2")
ax2.set_xlabel("Day of the Year")
ax2.set_ylim([0,2000])
ax2.set_title("Diffused")
ax2.set_xticks(x)
ax2.set_xticklabels(my_xticks)
ax2.plot(x,day_points[:,2],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax7=plt.twinx(ax2)
ax7.set_ylim([0,2000])
ax7.plot(x,day_points[:,5],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax8=plt.twinx(ax2)
ax8.set_ylim([0,2000])
ax8.plot(x,day_points[:,8],color='g',label="Emissivity vs. Wavelength",linewidth=2)

ax3=fig.add_subplot(313)
ax3.set_ylabel("Wh/m2")
ax3.set_xlabel("Day of the Year")
ax3.set_ylim([0,12000])
ax3.set_title("Global")
ax3.set_xticks(x)
ax3.set_xticklabels(my_xticks)
ax3.plot(x,day_points[:,3],color='r',label="Emissivity vs. Wavelength",linewidth=2)
ax9=plt.twinx(ax3)
ax9.set_ylim([0,12000])
ax9.plot(x,day_points[:,6],color='b',label="Emissivity vs. Wavelength",linewidth=2)
ax10=plt.twinx(ax3)
ax10.set_ylim([0,12000])
ax10.plot(x,day_points[:,9],color='g',label="Emissivity vs. Wavelength",linewidth=2)

fig.subplots_adjust(wspace=0.6, hspace=0.7)
fig.savefig("Solar_Radiation_Model_Comparision.png")


    
