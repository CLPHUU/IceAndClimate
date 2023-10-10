#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:37:18 2023

@author: berg
"""

import datetime as dt
import numpy as np
import SEB_functions as SEBf
import matplotlib.pyplot as plt

station = "S10"

if station == "S5":
    SEBdata = SEBf.SEB_data(FileName="PKM-data/S5_SEB_2003_2019_rp10b.txt")
    
elif station == "S6":
    SEBdata = SEBf.SEB_data(FileName="PKM-data/S6_SEB_2003_2019_rp4.txt")
    GsOld   = SEBdata.Extract_Variable("Gs")
    SEBdata.Correct_Gs_S6()

elif station == "S9":
    SEBdata = SEBf.SEB_data(FileName="PKM-data/S9_SEB_2003_2019_5.txt")
    
elif station == "S10":
    SEBdata = SEBf.SEB_data(FileName="PKM-data/S10_SEB_2009_2019.txt")


#%%
SWdown = SEBdata.Extract_Variable("SWin_corr")
SWup   = SEBdata.Extract_Variable("SWout")
SWnet  = SEBdata.Extract_Variable("SWnet_corr")
LWdown = SEBdata.Extract_Variable("LWin")
LWup   = SEBdata.Extract_Variable("LWout_corr")
LWnet  = SEBdata.Extract_Variable("LWnet_model")
SHF    = SEBdata.Extract_Variable("Hsen")
LHF    = SEBdata.Extract_Variable("Hlat")
Gs     = SEBdata.Extract_Variable("Gs")
MeltS  = SEBdata.Extract_Variable("melt_energy") 
MeltT  = SEBdata.Extract_Variable("totm_nrg") 
Residu = SEBdata.Extract_Variable("rest_energy")

if station == "S6":
    # recover the "apparrent GsRec
    SWint  = SEBdata.Extract_Variable("SumDivQ")
    GsRec  = Gs - SWint - MeltS + MeltT




#%%
# Example plots of the SEB
# OPTIONS HERE ARE
# Typical yearly cycle -> AvgMonth
# Monthly averages     -> Monthly
# Daily averages       -> Daily

PlotType = "AvgMonth"

if PlotType == "AvgMonth":
    # Example of plotting the typical yearly cycle
    MyFunc = SEBf.get_avg_monthly_value
    Range  = [0, 12.5]
    Label  = "Month"
    Xdata  = np.arange(14)
elif PlotType == "Monthly":
    # Example of montly averages
    MyFunc = SEBf.get_monthly_average
    Range  = [SEBdata.DateTime[0], SEBdata.DateTime[-1]]
    Label  = "Year"
    # do dummy call to function to get the date-times for the x-axis
    Xdata  = MyFunc(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]
elif PlotType == "Daily":
    # Example of montly averages
    MyFunc = SEBf.get_daily_average
    # restrain ourselves to one year, I take here 2016
    Range  = [dt.datetime.fromisoformat("2016-01-01"), dt.datetime.fromisoformat("2017-01-01")]
    Label  = "Date"
    # do dummy call to function to get the date-times for the x-axis
    Xdata  = MyFunc(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]


MonthSWdownS = MyFunc(SWdown, SEBdata.DateTime, PrintInfo=True)
MonthSWupS   = MyFunc(SWup, SEBdata.DateTime)

# however, it is easier to do the call to the SEBf-function in the plot call...

fig, axs = plt.subplots(2, sharex=True)

# upper figure are the radiative fluxes
axs[0].plot(Xdata, MonthSWdownS, 'b', linewidth=0.5, label="$SW_{down}$")
axs[0].plot(Xdata, -MonthSWupS, 'b:', linewidth=0.5, label="$SW_{up}$")
if station == "S6":
    axs[0].plot(Xdata, MyFunc(-SWint, SEBdata.DateTime), 'b--', linewidth=0.5, label="$SW_{int}$")

axs[0].plot(Xdata, MyFunc(LWdown, SEBdata.DateTime), 'r', linewidth=0.5, label="$LW_{down}$")
axs[0].plot(Xdata, MyFunc(-LWup, SEBdata.DateTime), 'r:', linewidth=0.5, label="$LW_{up}$")
axs[0].plot(Xdata, MyFunc(SWnet+LWnet, SEBdata.DateTime), 'k', linewidth=0.5, label="$R_{net}$")
if station == "S6":
    axs[0].plot(Xdata, MyFunc(SWnet+LWnet-SWint, SEBdata.DateTime), 'k:', linewidth=0.5, label="$R_{net surf}$")
axs[0].set_ylabel("Energy flux (W/m2)")
axs[0].legend(loc='lower right') # well, no spot is nice
axs[0].grid(True)

if station == "S6":
    axs[1].plot(Xdata, MyFunc(SWnet-SWint, SEBdata.DateTime), 'b', linewidth=0.5, label="$SW_{net surf}$")
else:
    axs[1].plot(Xdata, MyFunc(SWnet, SEBdata.DateTime), 'b', linewidth=0.5, label="$SW_{net}$")
axs[1].plot(Xdata, MyFunc(LWnet, SEBdata.DateTime), 'r', linewidth=0.5, label="$LW_{net}$")
axs[1].plot(Xdata, MyFunc(SHF, SEBdata.DateTime), 'seagreen', linewidth=0.5, label="$SHF$")
axs[1].plot(Xdata, MyFunc(LHF, SEBdata.DateTime), 'orange', linewidth=0.5, label="$LHF$")
axs[1].plot(Xdata, MyFunc(Gs, SEBdata.DateTime), 'grey', linewidth=0.5, label="$Gs$")
axs[1].plot(Xdata, MyFunc(MeltS, SEBdata.DateTime), 'purple', linewidth=0.5, label="$M$")
axs[1].plot(Xdata, MyFunc(Residu, SEBdata.DateTime), 'k', linewidth=0.5, label="Residue SEB model")

axs[1].set_ylabel("Energy flux (W/m2)")
axs[1].legend(loc='lower right') # Again, no spot is nice
axs[1].set_xlabel(Label)
axs[1].set_xlim(Range)
axs[1].grid(True)



#%%

# Example, a simple adjustment: Increase LW down by making the apparent atmospheric temperature 1 K warmer

dtemp = 1.
Tsurf = SEBdata.Extract_Variable("Tsurf_calc")
LWmod = SEBf.convert_T_in_LWout(Tsurf)

LWout = LWnet - LWdown
Tatm   = SEBf.convert_LWout_in_T(-LWdown, Celcius=False)
LWinDT = -SEBf.convert_T_in_LWout(Tatm+dtemp, Celcius=False)

LWoutDTuc = LWnet - LWinDT
TskinDT = SEBf.convert_LWout_in_T(LWoutDTuc)
TskinDT = np.where(TskinDT>0., 0., TskinDT)
LWoutDT = SEBf.convert_T_in_LWout(TskinDT)
MeltTDT = MeltT + LWoutDT - LWoutDTuc

fig2, ax2 = plt.subplots()
ax2.plot(SEBdata.DateTime, SEBf.get_running_melt_sum(MeltTDT, SEBdata.TimeStep), 'r', label='Melt for +1K')
ax2.plot(SEBdata.DateTime, SEBf.get_running_melt_sum(MeltT,   SEBdata.TimeStep), 'g', label='Observed Melt')
ax2.legend(loc='upper left')
ax2.set_ylabel('Accumulated melt (m w.e.)')
ax2.set_xlabel('Year')
ax2.set_xlim([SEBdata.DateTime[0], SEBdata.DateTime[-1]])
ax2.set_ylim([0, 19])




