#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 17:41:56 2023

@author: berg
"""

import datetime as dt
import numpy as np
import re
import os

# this file contains
# Methods to read and process data presented in 

# Authors: Kuipers Munneke, P., et al, 
#   Title: The K-transect on the western Greenland IceSheet: Surface energy balance (2003–2016)
# Journal: Arctic, Antarctic, and Alpine Research
#     DOI: 10.1080/15230430.2017.1420952

# the data themselves are received form P. Kuipers Munneke by personal communication.

# the included functions are
# the class SEB_data
# more info? -> help(SEB_data)


version = "1.1"
class SEB_data:
    '''Version 1.1 of a python class that reads and organizes SEB model output,
    derived from observations from the K-transect, Greenland'''
    

    def __init__(self, FileName=""):
        '''For the initialisation of this class, only the filename (including path) is needed.
        During the initialisation all data is read.
        The function works on all four provided data sets. It rectifies errors in the S5 and S6 data files.
        
        Output is a SEB_data-class object, containing the variables:
            ok             bool     data properly readed
            version        string   version of SEB_data
            FileName       string   name of the SEB data file, without path
            nval           integer  number of 
        '''
        
        self.ok = False
        self.version = "SEB_hourly_data version "+version
        self.FileName = re.split(r'/+', FileName)[-1] # remove path
        
        if FileName=="":
            print("A FileName is required!")
            return
        if not os.path.isfile(FileName):
            print("SEB data file "+FileName+" does not exist, return")
            return
        
        # first count number of lines and read the header
        with open(FileName) as AWSfile:
            self.nval = -1
            for line in AWSfile:
                if self.nval == -1:
                    # remove trailing line break and split line on the spaces
                    splitted_header = (line.strip()).split() 
                    print("The header has {0:3d} entries.".format(np.size(splitted_header)))

                    if splitted_header[0]=="year" and splitted_header[1]=="day":
                        if splitted_header[2]=="hour":
                            self.TimeStep = 3600.
                            self.Hourly   = True
                            ivstart = 3
                            ntime   = 3
                        else:    # assume daily data
                            self.TimeStep = 86400.
                            self.Hourly   = False
                            ivstart = 2
                            ntime   = 2
                    else:
                        print("We assumed that every header starts with 'year hour' (and a lot of spaces)")
                        print("However, this file starts with '"+splitted_header[0]+
                              "' and '"+splitted_header[1]+"'.")
                        print("Fix this! (or ask someone to fix this)")
                        return

                    if splitted_header[ivstart]=="Time":
                        ivstart+=1 # neglect this entry
                        
                    self.nvar = np.size(splitted_header) - ivstart
                    self.Variables = splitted_header[ivstart:]

                self.nval += 1
        print("AWS file '{0:s}' has {1:5d} lines of data for {2:2d} variables, start reading it.".format(self.FileName, 
                        self.nval, self.nvar))

        
        self.AllData    = np.zeros( [ self.nvar, self.nval ])                
        self.yyddhh     = np.zeros( [ 3, self.nval ], dtype=int )
        self.DateTime   = np.zeros( self.nval, dtype=dt.datetime )
        self.nvarDate   = np.zeros( self.nval, dtype=int)    # the number of variables per time entry. Should be constant, but isn't for AWS5 data :-(
        self.varvalid   = np.ones( self.nval, dtype=bool )
        
        # read data
        with open(FileName) as AWSfile:
            ival = -1
            for line in AWSfile:
                if ival>=0:
                    splitted_line = (line.strip()).split()
                    
                    # the time stap has errors, neglect is as a whole except first entry
                    self.yyddhh[:ntime,ival]   = [float(value) for value in splitted_line[:ntime]]

                    if ival==0:
                        # It is a bit lengthy to get the datetime format filled. 
                        # It does not help that timedelta needs int32, while the array is int64.
                        self.DateTime[ival]   = dt.datetime.fromisoformat(splitted_line[0][:4]+"-01-01") + dt.timedelta(days=int(self.yyddhh[1,ival]-1), hours=int(self.yyddhh[2,ival]))
                    else:
                        self.DateTime[ival] = self.DateTime[ival-1] + dt.timedelta(seconds=self.TimeStep)
                        
                    if self.yyddhh[0,ival] == -999:
                        lvalid = False
                        # fill year, doy and hour entries
                        self.yyddhh[0,ival] = self.DateTime[ival].year
                        ndays_dttype = self.DateTime[ival]-dt.datetime.fromisoformat(str(self.yyddhh[0,ival])+"-01-01")
                        self.yyddhh[1,ival] = ndays_dttype.days + 1
                        self.yyddhh[2,ival] = self.DateTime[ival].hour
                    else:
                        lvalid = True
                    
                    if lvalid==False:
                        self.AllData[:, ival] = np.zeros(self.nvar)*np.NAN
                        
                    elif np.size(splitted_line)-ivstart == self.nvar:
                        self.AllData[:, ival] = [float(value) for value in splitted_line[ivstart:]]
                        self.nvarDate[ival]   = self.nvar
                         
                    elif np.size(splitted_line)-ivstart == self.nvar+1 and self.FileName=="S6_SEB_2003_2019_rp4.txt":
# this needs to be fixed - it isn't fixable here. Specific data"repair" is done below.
                        self.AllData[:, ival] = [float(value) for value in splitted_line[ivstart:-1]]
                        self.nvarDate[ival]   = self.nvar+1

                    elif self.FileName=="S6_SEB_2003_2019_rp4.txt":
# this needs to be fixed - it isn't fixable here. Specific data"repair" is done below.
                        self.nvarDate[ival]   = np.size(splitted_line)-ivstart
                        self.AllData[:self.nvarDate[ival], ival] = [float(value) for value in splitted_line[ivstart:]]

                    else:
                        print("Errors should not occur in other dataset than the one for S6.")
                ival += 1

        
        print("Reading completed.")
        
        # turn invalid data into NaN
        self.AllData = np.where(self.AllData==-999., np.NaN, self.AllData)
        
        # Well, I wished the S6 data was error free and consistent, but it isn't.
        if self.FileName=="S6_SEB_2003_2019_rp4.txt":
            print("Apply data corrections...")
            # missed invalid data, likely due to station mainenance
            self.AllData[:,17270] = np.zeros(self.nvar)*np.NaN
            self.AllData[:,17271] = np.zeros(self.nvar)*np.NaN
            self.AllData[:,17273] = np.zeros(self.nvar)*np.NaN
            self.AllData[:,17274] = np.zeros(self.nvar)*np.NaN
            
            # a variable dissapears somehow...
            self.AllData[29:,114288:] = self.AllData[28:-1,114288:]
            self.AllData[28, 114288:] = self.AllData[28, 114288:] * np.NaN
            
            # variables 51 and higher look odd - remove
            self.AllData[51:,:] = self.AllData[51:,:] * np.NaN
            self.varvalid[51:]  = False
            print("We discard variable 51 ({0:s}, countinf from 0) and later as the data does not look consistent.".format(self.Variables[51]))
            print("It is further advised to call self.Correct_Gs_S6() as the Ground Heat Flux has sometimes erroneous data.")
        elif self.FileName=="S5_SEB_2003_2019_rp10b.txt" or self.FileName=="S9_SEB_2003_2019_5.txt" or self.FileName=="S10_SEB_2009_2019.txt":
            print("No data corrections needed for this station.")
        else:
            print("There are no data corrections know for this file, but that doesn't imply these are not necessary.")
                
        

        return
    
    def List_Variables(self):
        print("  #: Variable")
        for ivar in range(self.nvar):
            if self.varvalid[ivar]:
                print("{0:3d}:  {1:11s}".format(ivar+1, self.Variables[ivar]))
            else:
                print("{0:3d}: ERASED {1:11s}".format(ivar+1, self.Variables[ivar]))
    
    
    def Find_Variable(self, VarName):
        '''This function finds the entry number of a Variable in the AllData array. 
        Please note this function is case sensitive - your request should be exactly matching.'''
        
        lok      = True
        VarIndex = -1
        for v in range(self.nvar):
            if VarName.strip() == self.Variables[v].strip():
                VarIndex = v
        if VarIndex == -1:
            print("Variable '"+VarName+"' not found.")
            lok = False
        return VarIndex, lok  


    def Extract_Variable(self, VarName, NoDataForFail=True):
        '''This function extracts the data of a Variable from the main dataset AllData.
        Unless the optional parameter NoDataForFail is set to False, no data is 
          given if the requested variable name is not found.'''
        
        VarIndex, lok = self.Find_Variable(VarName)
        if not self.varvalid[VarIndex] and lok:
            print("WARNING: The data for variable '"+VarName+"' has been erased due to consistency problems.")
            print("         NO DATA IS PROVIDED")

        if lok and self.varvalid[VarIndex]:
            return self.AllData[VarIndex, :]
        elif NoDataForFail:
            return
        else:
            return np.zeros(self.nval)*np.NaN
        
    def Correct_Gs_S6(self):
        if self.FileName=="S6_SEB_2003_2019_rp4.txt":
            print("Recalculate the Ground Heat Flux by assuming that all other reported fluxes are right.")
            SWnet  = self.Extract_Variable("SWnet_corr")
            SWint  = self.Extract_Variable("SumDivQ")
            LWnet  = self.Extract_Variable("LWnet_model")
            SHF    = self.Extract_Variable("Hsen")
            LHF    = self.Extract_Variable("Hlat")
            MeltS  = self.Extract_Variable("melt_energy") 
            Residu = self.Extract_Variable("rest_energy")

            Gs = - SWnet + SWint - LWnet - SHF - LHF + Residu + MeltS 
            
            ivarGs = self.Find_Variable("Gs")
            self.AllData[ivarGs, :] = Gs

        else:
            print("This correction cannot be applied to other stations than S6")
        
        
def get_doy(year, month, day):
    timediff = dt.datetime(year=year, month=month, day=day) - dt.datetime(year=year, month=1, day=1)
    return timediff.days + 1

def convert_T_in_LWout(Tskin, Celcius=True):
    if Celcius:
        Tadd = 273.16
    else:
        Tadd = 0.
        
    return -( (Tskin+Tadd)**4 ) * 5.67E-8

def convert_LWout_in_T(LWout, Celcius=True):
    if Celcius:
        Tadd = 273.16
    else:
        Tadd = 0.
        
    return ( -LWout/5.67E-8 )**(0.25) - Tadd

def get_running_melt_sum(MeltE, TimeStep, ResetAtNan=True):
    '''The meltsum is in m w.e. per year'''
    RmeltSum = np.zeros(np.size(MeltE))
    for i in range(np.size(MeltE)):
        if not np.isnan(MeltE[i]):
            RmeltSum[i] = RmeltSum[i-1] + MeltE[i]*TimeStep/334000.
        elif not ResetAtNan:
            RmeltSum[i] = RmeltSum[i-1]
        elif np.isnan(MeltE[i-1]):
            RmeltSum[i-1] = np.NaN
    
    RmeltSum = RmeltSum/1000.
    
    return RmeltSum
            
def get_daily_average(Var, DateTime, GiveDatesBack=False, PrintInfo=False):
    '''This function derives the daily averages of a variable.
    It start at the first full day.'''
    Istart = (24-DateTime[0].hour)%24
    Ndays  = (DateTime[-1]-DateTime[0]).days
    
    if PrintInfo:
        print("The date and hour of the first entry is {}.".format(DateTime[0].strftime("%B %d, %Y, %H:%M")))
        print("The date and hour of the last entry is {}.".format(DateTime[-1].strftime("%B %d, %Y, %H:%M")))
        print("The dataset contains data of {0:5d} days.".format(Ndays))

# In order to do this in one call, I reshape the array into a 2D array and average over the second axis, thus the data of one day.    
    VarDay = np.mean(np.reshape(Var[Istart:Istart+Ndays*24], [Ndays, 24]),1)
    
    if GiveDatesBack:
        return VarDay, DateTime[Istart:Istart+Ndays*24:24]
    else:
        return VarDay
    
def get_next_month(DateTime):
    if DateTime.month==12:
        return dt.datetime(year=DateTime.year+1, month=1, day=1)
    else:
        return dt.datetime(year=DateTime.year, month=DateTime.month+1, day=1)
    
    
def get_monthly_average(Var, DateTime, GiveDatesBack=False, PrintInfo=False):
    '''This function derives the monthly averages of a variable.
    It start at the first full day.'''
    
    # find the start of the first full month
    NextMonth = get_next_month(DateTime[0])
    iStart = int( (NextMonth-DateTime[0]).days*24 + (NextMonth-DateTime[0]).seconds/3600 )

    # find the number of months
    nData  = np.size(DateTime)
    iMS    = iStart
    iME    = (get_next_month(DateTime[iMS])-DateTime[iMS]).days * 24
    nMonth = 0
    while iME<nData:
        iMS     = iME
        iME    += (get_next_month(DateTime[iME])-DateTime[iME]).days * 24
        nMonth += 1
    iFinal = iMS    
    
    if PrintInfo:
        print("The date and hour of the first entry is {}.".format(DateTime[0].strftime("%B %d, %Y, %H:%M")))
        print("First full month starts at {}".format(DateTime[iStart].strftime("%B %d, %Y, %H:%M")))
        print("The date and hour of the last entry is {}.".format(DateTime[-1].strftime("%B %d, %Y, %H:%M")))
        print("Last full month end at {}".format(DateTime[iFinal-1].strftime("%B %d, %Y, %H:%M")))
        print("The dataset contains data of {0:3d} months.".format(nMonth))
    
    
    VarMonth = np.zeros(nMonth)
    VarDate  = np.zeros(nMonth, dtype=dt.datetime )
    iMS = iStart
    for iMonth in range(nMonth):
        iME = iMS + (get_next_month(DateTime[iMS])-DateTime[iMS]).days * 24
        VarMonth[iMonth] = np.mean(Var[iMS:iME])
        VarDate[iMonth]  = DateTime[iMS]
        iMS = iME
        
    if GiveDatesBack:
        return VarMonth, VarDate
    else:
        return VarMonth
    
    
def get_avg_monthly_value(Var, DateTime, PrintInfo=False, GiveNumberOfMonths=False):
    '''This function derives the mean monthly value of a variable'''
    VarMonth, VarDate = get_monthly_average(Var, DateTime, GiveDatesBack=True, PrintInfo=PrintInfo)
    
    VarMeanMonth = np.zeros(14)
    Ndata        = np.zeros(13)
    for iMonth in range(np.size(VarMonth)):
        iM = VarDate[iMonth].month
        # leave out months with no valid data
        if not np.isnan(VarMonth[iMonth]):
            VarMeanMonth[iM] += VarMonth[iMonth]
            Ndata[iM]        += 1
    
    VarMeanMonth[1:13] = VarMeanMonth[1:13]/Ndata[1:13]
    VarMeanMonth[0]    = VarMeanMonth[12]
    VarMeanMonth[13]   = VarMeanMonth[1]
    
    if GiveNumberOfMonths:
        return VarMeanMonth, Ndata
    else:
        return VarMeanMonth
   


