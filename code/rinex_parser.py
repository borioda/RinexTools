#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 17:20:26 2022

@author: daniele
"""

import numpy as np
from tqdm.auto import tqdm
import os

import timefun as tf
import gnss_ephemeris as ge
from datetime import datetime


class nav_rinex_parser :
    
    # Add the codes corresponding to the ephemerides implemented
    gnss_codes = ['G', 'E']
    
    
    """
    Summmary :
        Class that parses a rinex 3 navigation file
    """
    def __init__(self, filename : str) :
        """
        Object constructor.

        Parameters
        ----------
        filename : str
            Path of the rinex navigation file

        Returns
        -------
            The nav rinex parser
        """
        self.filename = filename
        
        # Try to open the file
        self.rinex_file = open(filename, 'r')
        
        # Now read parse the header
        self.header = self.parse_header()
        
    
    def parse_header(self) :
        """
        Extract the lines corresponding to the headers. 
        For the moment these lines are stored in a list.

        Returns
        -------
        List of lines corresponding to the header.
        """
        header = []
        
        line = self.rinex_file.readline()
        
        while "END OF HEADER" not in line :
            header.append(line)
            line = self.rinex_file.readline()
    
    def get_ephemeris(self) :
        
        # Read lines until a valid code is found 
        line = self.rinex_file.readline()
        
        if "" == line :
            return None
        
        while line[0] not in nav_rinex_parser.gnss_codes :
            line = self.rinex_file.readline()
            
            if "" == line :
                return None
            
        # If here a valid GNSS code has been found
        gnss = line[0]
        if gnss == 'G' or gnss == 'E' :
            
            # Build a data block of 8 lines
            data_block = [line]
            for ii in range(7) :
                line = self.rinex_file.readline()
                
                if "" == line :
                    return None
                
                data_block.append(line)
                
        # Create the corresponding ephemeris
        eph = None
        
        if gnss == 'G' :
            eph = ge.gps_rinex_ephemeris(data_block)
            
        if gnss == 'E' :
            eph = ge.gal_rinex_ephemeris(data_block)
            
        return eph
    
class rinex_parser :
    """
    Summary:
        Class that reads a measurement rinex file and convert it into a csv
        file
    """
    
    def __init__(self, filename : str) :
        """
        Summary:
            Object constructor.
        """ 
        
        self.rinex_file_name = filename
        
        # open the rinex file
        self.rinex_file = open(filename, 'r')
        
        # read the header
        header = []
        line = self.rinex_file.readline()
        
        while "END OF HEADER" not in line :
            header.append(line)
            line = self.rinex_file.readline()
            
        # interpret the header
        self.obs_list = rinex_parser.get_obs_list(header)
        
        self.first_obs_time = rinex_parser.get_time_first_obs(header)
        
        self.last_obs_time = rinex_parser.get_time_last_obs(header)
        
        self.obs_interval = rinex_parser.get_obs_interval(header)
    
    @staticmethod
    def get_obs_list(header : list) -> dict :
        """
        Summary:
            Extract the list of observations from the rinex header.
        
        Arguments:
            header : list of lines representing the header.
            
        Return :
            Dictionary with the list of observations associated to the different
            GNSSs.
        """
        return_dict = {}
        gnss = 'G'
        
        for line in header :
            if "SYS / # / OBS TYPES" in line :
                if line[0] in "GECRS" :
                    gnss = line[0]
                    return_dict[gnss] = line[:60].split()[2:]
                else :
                    return_dict[gnss] = [*return_dict[gnss], *line[:60].split()] 
        
        return return_dict
    
    @staticmethod
    def get_time_first_obs(header : list) -> datetime :
        
        """
        Summary:
            Extract the datetime of the first observation is present in the header.
        
        Arguments:
            header : list of lines representing the header.
            
        Return :
            datetime of the first observation.
        """
        starttime = None
        
        for line in header :
            if "TIME OF FIRST OBS" in line :
                dateinfo = line.split()
                
                seconds = float(dateinfo[5]) 
                microsec = int((seconds - np.floor(seconds)) * 1e6)
                seconds = int(np.floor(seconds))
                
                starttime = datetime(int(dateinfo[0]), int(dateinfo[1]), int(dateinfo[2]), \
                                     hour = int(dateinfo[3]), minute = int(dateinfo[4]), \
                                     second = seconds, microsecond = microsec )
                break
        
        return starttime
    
    @staticmethod
    def get_time_last_obs(header : list) -> datetime :
        
        """
        Summary:
            Extract the datetime of the last observation is present in the header.
        
        Arguments:
            header : list of lines representing the header.
            
        Return :
            datetime of the last observation.
        """
        stoptime = None
        
        for line in header :
            if "TIME OF LAST OBS" in line :
                dateinfo = line.split()
                
                seconds = float(dateinfo[5]) 
                microsec = int((seconds - np.floor(seconds)) * 1e6)
                seconds = int(np.floor(seconds))


                stoptime = datetime(int(dateinfo[0]), int(dateinfo[1]), int(dateinfo[2]), \
                                     hour = int(dateinfo[3]), minute = int(dateinfo[4]), \
                                     second = seconds, microsecond = microsec )
                break
        
        return stoptime
    
    @staticmethod
    def get_obs_interval(header : list) -> float :    
        """
        Summary:
            Extract observation interval from the rinex header.
        
        Arguments:
            header : list of lines representing the header.
            
        Return :
            The observation interval in seconds
        """
        obsint = 1.0
        
        for line in header :
            if "INTERVAL" in line :
                obsint = float(line.split()[0])
                break
        
        return obsint
    
    def get_obsnumber(self) :
        
        if (self.first_obs_time is not None) and (self.last_obs_time is not None) :
            obsnum = (self.last_obs_time - self.first_obs_time).total_seconds() / self.obs_interval
            
            return obsnum 
    
    def to_csv(self, filename : str) :
        """
        Summary:
            Convert the rinex file to a csv file
        
        Arguments:
            filename - output csv file.
            
        Return :
            Nothing.
        """
        output_file = open(filename, 'w')
        output_file.write("TOW,WEEK,GNSS,PRN,OBS_TYPE,VALUE\n")
        
        with tqdm(total=os.path.getsize(self.rinex_file_name)) as pbar:
            # get current position in file
            pbar.update(self.rinex_file.tell())
            
            while True :
                line = self.rinex_file.readline()
                if "" == line :
                    break
                
                # This line should contain the date and the number of observations
                tow, week, nobs = rinex_parser.get_obs_info(line)
                
                data_block = []
                for ii in range(nobs) :
                    line = self.rinex_file.readline()
                    if "" != line :
                        data_block.append(line)
                
                self.print_obs_block(output_file, tow, week, data_block)
                
                pbar.update(self.rinex_file.tell() - pbar.n)
                    
        output_file.close()
                
    def print_obs_block(self, output_file, tow : float, week : int, data_block : list) :
        
        for line in data_block :
            gnss = line[0]
            prn = int(line[1:3])
            
            # get the list of observations
            if gnss in self.obs_list :
                obs_list = self.obs_list[gnss]
            else :
                # observation not supported
                return
            
            # take into account empty spaces
            obs = line.replace(' '*16, ' -1').split()[1:]
            
            # remove signal strength indicators
            obs_final = [x for x in obs if len(x) > 1] 
            
            if len(obs_final) > len(obs_list) :
                return
            
            # finally print the observations
            for ii, obs in enumerate(obs_final) :
                if obs != "-1" :
                    output_file.write(f'{tow},{week},{gnss},{prn},{obs_list[ii]},{obs}\n')
    
    def get_obs_block(self, offset = 0) :
        
        # Read lines until a symbol '>' is found
        line = self.rinex_file.readline()
        
        if "" == line :
            return None
        
        while line[0] != '>' :
            line = self.rinex_file.readline()
            
            if "" == line :
                return None
            
        # Now line contains the epoch information
        
        # Build the observation block
        info = line.split()
        nobs = int(info[-1 - offset])
        obs_block = {'year' : int(info[1]),
                     'month' : int(info[2]),
                     'day' : int(info[3]),
                     'hour' : int(info[4]),
                     'minutes' : int(info[5]),
                     'seconds' : float(info[6]),
                     'obs' : [] 
                     }
        
        # Now get the actual observations
        for ii in range(nobs) :
            line = self.rinex_file.readline()
            
            gnss = line[0]
        
            # get the list of observations
            if gnss in self.obs_list :
                obs_list = self.obs_list[gnss]
            else :
                # observation not supported
                continue
            
            obs = rinex_parser.parse_obs(line, obs_list)
            
            if obs is not None :
                obs_block['obs'].append(obs.copy())
            
        return obs_block
    

            
    @staticmethod
    def parse_obs(line : str, obs_list) :
        
        gnss = line[0]
        prn = int(line[1:3])
        
        # take into account empty spaces
        obs = line.replace(' '*16, ' -1').split()[1:]
        
        # remove signal strength indicators
        obs_final = [x for x in obs if len(x) > 1] 
        
        # There may be missing observations at the end!
        # So, obs_final cannot be longer than obs_list
        if len(obs_final) > len(obs_list) :
            return None
        
        obs_dict = {
            'gnss' : gnss,
            'prn' : prn
            }
        
        for ii, val in enumerate(obs_final) :
            
            key = obs_list[ii]
            
            if val != "-1" :
                obs_dict[key] = float(obs_final[ii])
        
        return obs_dict
        
    @staticmethod
    def get_obs_info(line : str) :
        info = line.split()
        
        nobs = int(info[-1])
        
        tow, week = tf.DateToGPS( int(info[1]), int(info[2]), int(info[3]), int(info[4]) )
        tow += int(info[5]) * 60 + float(info[6])
        
        return tow, week, nobs


def print_block(obs_dict : dict, rinex_file, obs_block : dict) :
    """
    Print an observation block to rinex file.
    rinex_file is an file pointer opened in write mode.
    """
    
    # First print the date
    line = "> " + str(obs_block["year"]) + " "
    if obs_block["month"] < 10 :
        line = line + " 0" + str(obs_block["month"])
    else :
        line = line + " " + str(obs_block["month"])
        
    if obs_block["day"] < 10 :
        line = line + " 0" + str(obs_block["day"])
    else :
        line = line + " " + str(obs_block["day"])
    
    if obs_block["hour"] < 10 :
        line = line + " 0" + str(obs_block["hour"])
    else :
        line = line + " " + str(obs_block["hour"])
        
    if obs_block["minutes"] < 10 :
        line = line + " 0" + str(obs_block["minutes"])
    else :
        line = line + " " + str(obs_block["minutes"])
        
    if obs_block["seconds"] < 10 :
        line = line + "  " + f'{obs_block["seconds"]:.7f}  0 '
    else :
        line = line + " " + f'{obs_block["seconds"]:.7f}  0 '
        
    # Now add the number of observations
    nobs = len(obs_block["obs"])
    
    if nobs < 10 :
        line = line + " " + str(nobs) + "\n"
    else :
        line = line + str(nobs) + "\n"
        
    # Write info to fine
    rinex_file.write(line)
    
    # Now print the actual measurements
    for obs in obs_block["obs"] :
        gnss = obs["gnss"]
        
        # list with the order of the observations
        if gnss in obs_dict :
            obs_list = obs_dict[gnss]
        else :
            continue
        
        # Add gnss and prn
        line = obs["gnss"] 
        if obs["prn"] < 10 :
            line = line + "0" + str(obs["prn"]) + " "
        else :
            line = line + str(obs["prn"]) + " "
        
        for key in obs_list :
            if key in obs :
                val_as_str = f'{obs[key]:.5f}'
            else :
                val_as_str = ''
                
            # pad with blanks
            if len(val_as_str) < 16 :
                val_as_str = ' '*(16 - len(val_as_str)) + val_as_str
                
            line = line + val_as_str
                
        line = line + "\n"
        
        rinex_file.write(line)

def add_meas_header( obs_dict : dict, rinex_file ):
    
    rinex_file.write('     3.02           OBSERVATION DATA    M: MIXED            RINEX VERSION / TYPE\n')
    rinex_file.write(' '.rjust(60) + 'MARKER NAME\n')
    rinex_file.write(' '.rjust(60) + 'MARKER NUMBER\n')
    rinex_file.write(' '.rjust(60) + 'MARKER TYPE\n')
    rinex_file.write(' '.rjust(60) + 'OBSERVER / AGENCY\n')
    rinex_file.write(' '.rjust(60) + 'REC # / TYPE / VERS\n')
    rinex_file.write(' '.rjust(60) + 'ANT # / TYPE\n')
    rinex_file.write('        0.0000        0.0000        0.0000'.ljust(60) + 'APPROX POSITION XYZ\n')
    rinex_file.write('        0.0000        0.0000        0.0000'.ljust(60) + 'ANTENNA: DELTA H/E/N\n')

    for gnss, obs_list in obs_dict.items() :
        line = gnss + '    '
        nobs = len(obs_list)
        
        if nobs < 10 :
            line = line + ' ' + str(nobs) + ' '
        else :
            line = line + str(nobs) + ' '
        
        for ii in range(min([13, nobs])) :
            line = line + obs_list[ii] + ' '
        
        line = line + ' '*( 60 - len(line)) + 'SYS / # / OBS TYPES\n' 
        
        # now write to file
        rinex_file.write(line)
        
        if nobs > 13 :
            line = ' '*7
            
            for ii in range(13, nobs) :
                line = line + obs_list[ii] + ' '
                
            # pad with white spaces
            line = line + ' '*( 60 - len(line)) + 'SYS / # / OBS TYPES\n'
            
            rinex_file.write(line)
    
    rinex_file.write(' '.rjust(60) + 'END OF HEADER\n')