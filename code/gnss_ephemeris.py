#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 09:19:58 2021

@author: daniele
"""

import abc
import numpy as np
import timefun as tf

class base_ephemeris(metaclass = abc.ABCMeta) :
    
    MU_GPS = 3.9860050e14     
    MU_GLO = 3.9860044e14
    MU_GAL = 3.986004418e14
    
    OMEGA_E = 7.2921151467e-5
    """
    Summary:
    
    """
    def __init__(self, gnss : str, prn : np.uint8) -> None : 
        """ 
        Summary:
            Object constructor.
            
        Arguments:
            gnss - the GNSS identifier
            prn - the satellite identifier (prn)
        Returns: 
            Nothing.
        """
        self.prn = prn
        self.gnss = gnss
      
    @abc.abstractmethod
    def get_sat_pos_and_vel( self, time : float ) :
        """
        Summary: 
            Get the satellite position and velocity at a given time.
            Satellite position and velocity are expressed in ECEF

        Arguments:
            time - transmission time for the position/velocity computation

    	Returns:
            Array with the satellite position and velocity    
        """
        pass
    
    @abc.abstractmethod
    def get_clock_corrections( self, time : float ) :
        """ 
        Summary:
            Get the satellite clock corrections including relativistic effects.
	
        Arguments:
            time - the time of computation
	
        Returns:
            The satellite clock bias and clock drift.
        """
        pass
    
    @abc.abstractmethod 
    def get_refence_time( self ) :
        pass
    
    @staticmethod
    def time_diff( reftime : float, time : float ) -> float :
        """ 
        Summary:
            Compute the absolute time difference between reftime and time.
            Possible week cross-overs are accounted for.
        
        Arguments:
            reftime - reference time
            time - the time with respect to which time differences are computed

        Returns:
            The absolute time difference
        """
        time_diff = time - reftime
        
		# check for end of week crossover
        if time_diff < -302400.0 :
            time_diff += 604800.0
		
        if time_diff > 302400.0 :
            time_diff -= 604800.0
        
        return time_diff
    
    def __eq__(self, other) :
        
        is_equal = isinstance(other, base_ephemeris)
        
        # Now check other conditions
        is_equal = is_equal and (self.prn == other.prn)
        is_equal = is_equal and (self.gnss == other.gnss)
        is_equal = is_equal and (self.get_refence_time() == other.get_refence_time())
        
        return is_equal
    
    def __str__(self):
        return self.__dict__.__str__()


class sp3_ephemeris(base_ephemeris) :
    
    """
    Summary :
        Class implementing precise ephemeris loaded from SP3 files
    """
    def __init__(self, gnss : str, prn : int, eph_info : dict ) :
        """
        Summary:
            Object constructor.
            
        Arguments:
            gnss - the GNSS identifier
            prn - the satellite identifier (prn)
            eph_info - dictionary with additional information regarding the
                       ephemeris
                       
                      it has the following fields:
                          eph_info = {'first_epoch_tow' : ,
                                      'week_number' :,
                                      'num_of_epochs' :,
                                      'epoch_interval' :,
                                      'coord_sys':,
                                      'accuracy':,
                                      'base_for_orb':,
                                      'base_for_clk':
                                    }
        """
        # base ephemeris
        super().__init__(gnss, prn)
        
        # save the ephemeris information
        self.eph_info = eph_info
        
        # allocate the vector with time and corrections
        self.times_of_corr = np.arange(eph_info['first_epoch_tow'], \
                                       eph_info['first_epoch_tow'] + \
                                       eph_info['num_of_epochs'] * eph_info['epoch_interval'], \
                                       eph_info['epoch_interval'] )
            
        # allocate the matrix with the corrections
        self.corr = np.zeros((len(self.times_of_corr), 4), dtype = float)
    
    def add_correction(self, time_of_corr, sp3_line : str) :
        """
        Summary :
            Function that add a correction to the correction array from a line
            of an sp3 file.
        
        Arguments :
            time_of_corr : the time of the correction. Either a single tow or
                           a dictionary with the explicit time
            sp3_line : string with the sp3 line
        
        Returns:
            Nothing.
        """
        # Check if this the correct GNSS and PRN
        if sp3_line[0] != 'P' :
            return
        
        # Correct GNSS
        if sp3_line[1] != self.gnss :
            return
        
        # Correct PRN
        if int(sp3_line[2:4]) != self.prn :
            return
        
        # Determine the time of the correction
        if isinstance(time_of_corr, dict) :
            tow, gps_week = tf.DateToGPS( time_of_corr['year'], \
                                          time_of_corr['month'], \
                                          time_of_corr['day'], \
                                          time_of_corr['hour'] )
            tow += time_of_corr['minute'] * 60 + time_of_corr['second']
        else :
            tow = time_of_corr
            
        # Now check if the tow is in the times of corrections
        ind = np.argwhere(self.times_of_corr == tow).flatten()
        
        if len(ind) == 0 :
            return
        else :
            ind = ind[0]
        
        # if here we can populate the corrections
        self.corr[ind][0] = float(sp3_line[4:18])
        self.corr[ind][1] = float(sp3_line[18:32])
        self.corr[ind][2] = float(sp3_line[32:46])
        self.corr[ind][3] = float(sp3_line[47:60])
                    
    def get_clock_corrections( self, time : float ) :
        
        return 0
        
    def get_refence_time(self) :
        return 0
    
    def get_sat_pos_and_vel( self, time : float ) :
        pass


class base_rinex_ephemeris(base_ephemeris) :
        
    def __init__(self, data_block : list) :
        """
        Summary:
            Initialize a generic ephemeris to be parsed from a rinex file
        """
        
        # star with the first line
        line = data_block[ 0 ]
        
        # get the GNSS ID
        gnss = line[0]
            
        # get the prn
        prn = np.uint8(line[1:3])
        
        # call the base class
        super().__init__(gnss, prn)
        
        # Member objects
        self.toc = float(0)
        self.week_number = np.uint16(0) 
        
        self.a_f0 = float(0)
        self.a_f1 = float(0)
        self.a_f2 = float(0) 
        
        self.iode = np.int16(0)
        self.crs = float(0)
        self.deltan = float(0)
        self.m0 = float(0)   
        
        self.cuc = float(0)
        self.e = float(0)
        self.cus = float(0)        
        self.sqrtA = float(0)
        
        self.toe = float(0)
        self.cic = float(0)
        self.omega0 = float(0)
        self.cis = float(0)
        
        self.i0 = float(0)
        self.crc = float(0)
        self.omega = float(0)
        self.omegaDot = float(0)
        
        self.iDot = float(0)
        
        self.accuracy = float(0)
        self.health = np.uint8(0)
        self.tgd = float(0)
        self.iodc = np.uint16(0)
        
        self.tow = float(0)
        self.fitInterval = float(0)
        
        # Now parse the data block
        self.parse_rinex_block(data_block)
        
        self.MU = float(0)

    def parse_rinex_block( self, data_block : list) :
        
        # star with the first line
        line = data_block[ 0 ]
        
        # get the time of clock
        year = np.uint16( line[4:8] )
        month = np.uint8( line[9:11] )
        day = np.uint8( line[12:14] )
        hour = np.uint8( line[15:17] )
        minutes = np.uint8( line[18:20] )
        seconds = np.uint8( line[21:23] )
        
        # Clock reference time in seconds and week number
        self.toc, self.week_number = tf.DateToGPS(year, month, day, hour)
        self.toc += 60 * minutes + seconds
        
        # Now get the clock parameters
        self.a_f0 = float(line[23:42].replace('D', 'e'))
        self.a_f1 = float(line[42:61].replace('D', 'e'))
        self.a_f2 = float(line[61:80].replace('D', 'e')) 
        
        # second line
        line = data_block[1]
        self.iode = np.int16(float(line[4:23].replace("D","e")))
        self.crs = float(line[23:42].replace("D","e"))
        self.deltan = float(line[42:61].replace("D","e"))
        self.m0 = float(line[61:80].replace("D","e"))   
        
        # third line
        line = data_block[2]
        self.cuc = float(line[4:23].replace("D","e"))
        self.e = float(line[23:42].replace("D","e"))
        self.cus = float(line[42:61].replace("D","e"))
        self.sqrtA = float(line[61:80].replace("D","e"))
        
        # fourth line
        line = data_block[3]
        self.toe = float(line[4:23].replace("D","e"))
        self.cic = float(line[23:42].replace("D","e"))
        self.omega0 = float(line[42:61].replace("D","e"))
        self.cis = float(line[61:80].replace("D","e"))
        
        # fifth line
        line = data_block[4]
        self.i0 = float(line[4:23].replace("D","e"))
        self.crc = float(line[23:42].replace("D","e"))
        self.omega = float(line[42:61].replace("D","e"))
        self.omegaDot = float(line[61:80].replace("D","e"))
        
        # sixth line
        line = data_block[5]
        self.iDot = float(line[4:23].replace("D","e"))
        # code on L2, GPS week and L2 P flag not read (GPS week already computed)
        
        # seventh line
        line = data_block[6]
        self.accuracy = float(line[4:23].replace("D","e"))
        self.health = np.uint8(float(line[23:42].replace("D","e")))
        self.tgd = float(line[42:61].replace("D","e"))
        self.iodc = np.uint16(float(line[61:80].replace("D","e")))
        
        # Last line
        line = data_block[7]
        self.tow = float(line[4:23].replace("D","e"))
        
        if line[23:42].replace("D","e").isnumeric() :
            self.fitInterval = float(line[23:42].replace("D","e"))
      
    
    def get_sat_pos_and_vel( self, time : float ) :
        """
        Summary: 
            Get the satellite position and velocity at a given time.
            Satellite position and velocity are expressed in ECEF

        Arguments:
            time - transmission time for the position/velocity computation

    	Returns:
            Array with the satellite position and velocity    
        """
        A = self.sqrtA**2
        
        # First compute the delta time tk
        tk = base_ephemeris.time_diff( time, self.toe )
        
        # Mean motion
        n0 = np.sqrt( self.MU / A**3 )
        
        # Corrected mean motion
        n = n0 + self.deltan

        # Compute the mean anomaly for tk
        mk = self.m0 + n * tk
        mk = ( mk + 2*np.pi ) % (2*np.pi)
                    
        # Solve iteratively for the Kepler equation for the eccentric anomaly Ek
        ek = mk
        ek_dot = n 
        
        for ii in range(10) :
            ek = mk  + self.e * np.sin(ek)
            ek_dot = n + self.e * np.cos(ek) * ek_dot
              
        ek = (ek + 2*np.pi) % (2*np.pi)
        
        roote = np.sqrt(1 - self.e**2)
        
        # Compute the true anomaly, vk
        sinEk = np.sin(ek)
        cosEk = np.cos(ek)

        vk = np.arctan2( roote * sinEk, cosEk - self.e )
                
        sin2 = np.sin(2*(self.omega + vk))
        cos2 = np.cos(2*(self.omega + vk))
        
        # Compute the argument of latitude
        uk = vk + self.omega + self.cuc * cos2 + self.cus * sin2
                             
        # Compute the radial distance
        rk = A * (1 - self.e * cosEk) + self.crc * cos2 + self.crs * sin2
        
        # Compute the inclination
        ik = self.i0 + self.iDot * tk + self.cic * cos2 + self.cis * sin2
        
        # Compute the longitude of the ascending node
        lambdak = self.omega0 + (self.omegaDot - base_ephemeris.OMEGA_E) * tk \
                              - base_ephemeris.OMEGA_E * self.toe
        
        # Compute the satellite position
        pos = np.zeros(3)
        x = rk * np.cos(uk) 
        y = rk * np.sin(uk) 
        cosi = np.cos(ik)
        sini = np.sin(ik)
        
        cosl = np.cos(lambdak)
        sinl = np.sin(lambdak)
        
        pos[0] = x * cosl - y * cosi * sinl
        pos[1] = x * sinl + y * cosi * cosl
        pos[2] = y * sini

        # Compute the velocity
        vel = np.zeros(3)
        lambdak_dot = self.omegaDot - base_ephemeris.OMEGA_E
        
        uk_dot = n * roote * (A / rk)**2
        
         #....... 4.2  compute the derivative of the 2nd harmonic corrections
         # with respect to time tk. (i.e. compute the derivative of:
         # 1 - argument of latitude correction
         # 2 - radius correction
         # 3 - correction to inclination)

        sin2pd = 2.0 * cos2 * uk_dot
        cos2pd = -2.0 * sin2 * uk_dot
        corlatd = self.cuc * cos2pd + self.cus * sin2pd
        corrd = self.crc * cos2pd + self.crs * sin2pd
        corid = self.cic * cos2pd + self.cis * sin2pd
        
        
        #....... 4.3  compute the derivative of the following w.r.t. time tk:
        #  1) corrected argument of latitude
        #  2) corrected radius
        #  3) corrected inclination  

        uk_dot += corlatd
        rk_dot = A * self.e * sinEk * ek_dot + corrd
        xikd = corid        
        
        # in plane velocities
        x_dot = rk_dot * np.cos(uk) - y * uk_dot
        y_dot = rk_dot * np.sin(uk) + x * uk_dot
        
        vel[0] = - x * sinl * lambdak_dot - y * cosl * cosi * lambdak_dot \
                 + y * sini * xikd * sinl \
                 + x_dot * cosl - y_dot * cosi * sinl
        vel[1] = x * cosl * lambdak_dot - y * sinl * cosi * lambdak_dot \
                 - y * sini * xikd * cosl\
                 + x_dot * sinl + y_dot * cosi * cosl
        vel[2] = y_dot * sini + y * cosi * xikd
                 
        return pos, vel
    
    def get_clock_corrections( self, time : float ) :
        
        return 0
        
    def get_refence_time(self) :
        return self.toc

class gps_rinex_ephemeris(base_rinex_ephemeris) :
    
    def __init__(self, data_block : list) :
        """
        Summary:
            Initialize a GPS ephemeris using the 8 lines data block from a rinex
            file
            
        Arguments:
            data_block - block of 8 lines with the data from a rinex file
            
        Returns:
            Nothing.
        """
        
        # call the base class
        super().__init__(data_block)
        
        # Check if the GNSS is the correct one          
        if self.gnss != 'G' :
            raise Exception("gps_rinex_ephemeris - __init__() - Not a GPS data block")
                    
        self.MU = base_ephemeris.MU_GPS
        
class gal_rinex_ephemeris(base_rinex_ephemeris) :
    
    def __init__(self, data_block : list) :
        """
        Summary:
            Initialize a Galileo ephemeris using the 8 lines data block from a rinex
            file
            
        Arguments:
            data_block - block of 8 lines with the data from a rinex file
            
        Returns:
            Nothing.
        """
        
        # call the base class
        super().__init__(data_block)
        
        # Check if the GNSS is the correct one          
        if self.gnss != 'E' :
            raise Exception("galileo_rinex_ephemeris - __init__() - Not a Galileo data block")
        
        # Also read the BGDb (not available in GPS)
        line = data_block[6]
        self.bgdb = float(line[61:80].replace("D","e"))
        
        self.MU = base_ephemeris.MU_GAL
    
        
class ephemeris_factory :
    
    def create_from_sp3(self, filename : str) :
        # first open the sp3 file
        sp3_file = open(filename, 'r')
        
        # Now interpret the header
        # First line
        # Example: #dP2021  5 23  0  0  0.00000000      96 d+D   IGb14 FIT AIUB
        line = sp3_file.readline()
        if line[0] != '#' :
            return
        
        # For the moment only position information is supported
        if line[2] != 'P' :
            return
        
        # In this cae we do not use the time information which is obtained from
        # the second line in terms of GPS week and tow
        eph_info = {}
        eph_info['num_of_epochs'] = int(line[32:39])
        eph_info['coord_sys'] = line[46:51]
        
        # Now move to the second line
        line = sp3_file.readline()
        
        eph_info['week_number'] = int(line[3:7])
        eph_info['first_epoch_tow'] = float(line[8:23])
        eph_info['epoch_interval'] = float(line[24:38])
        
        # Now determine the satellites supported in the file
        line = sp3_file.readline()
        numSat = int(line[4:6])
        
        sat_list = [line[ii:(ii+3)] for ii in range(9,60,3)]
        for ii in range(4) :
            line = sp3_file.readline()
            if len(sat_list) < numSat :
                sat_list = [*sat_list, *[line[ii:(ii+3)] for ii in range(9,60,3)]]

        sat_list = sat_list[:numSat]            
        
        # Do the same for the satellite accuracy
        line = sp3_file.readline()
        
        acc_list = [int(line[ii:(ii+3)]) for ii in range(9,60,3)]
        for ii in range(4) :
            line = sp3_file.readline()
            if len(acc_list) < numSat :
                acc_list = [*acc_list, *[int(line[ii:(ii+3)]) for ii in range(9,60,3)]]

        acc_list = sat_list[:numSat]
        
        # skip two lines (it should be used to get the system time)
        line = sp3_file.readline()
        line = sp3_file.readline()
        
        # now get the base for Pos/Vel and Clk/rate
        line = sp3_file.readline()
        eph_info['base_for_orb'] = float(line[3:13])
        eph_info['base_for_clk'] = float(line[14:26])
        
        # skip other lines
        for ii in range(3) :
            line = sp3_file.readline()
            
        # Now skip comments
        line = sp3_file.readline()
        while line[0:2] == '/*' :
            line = sp3_file.readline()
            
        # if here the header is over :
        # create the return list for the ephemeris
        eph_list = []
        for ii in range(numSat) :
            eph_info['accuracy'] = acc_list[ii]
            new_eph = sp3_ephemeris(sat_list[ii][0], int(sat_list[ii][1:3]), eph_info )
            eph_list.append(new_eph)
            
        while line != 'EOF' :
            time = ephemeris_factory.get_time_from_sp3_line(line)
            
            if time is None :
                break
            
            for ii in range(numSat) :
                line = sp3_file.readline()
                eph_list[ii].add_correction(time, line)
            
            # eventally read the new date
            line = sp3_file.readline()
        
        return eph_list
        
    @staticmethod
    def get_time_from_sp3_line(sp3_line : str) -> dict :
        
        if sp3_line[0] != '*' :
            return None
        
        date = {'year' : int(sp3_line[3:7]),
                'month' : int(sp3_line[8:10]),
                'day' : int(sp3_line[11:12]),
                'hour' : int(sp3_line[14:16]),
                'minute' : int(sp3_line[17:19]),
                'second' : float(sp3_line[20:31])
            }
        
        return date
    
    def create_from_rinex( self, filename : str ) :
        
        # first open the rinex file
        rinex_file = open(filename, 'r')
        
        # create the output list
        eph_list = []
        
        # start reading the rinex file line by line
        lines = []
        line = rinex_file.readline()
        while 'END OF HEADER' not in line :
            lines.append(line)
            line = rinex_file.readline()
            
        # interpret the header
        # header_info = interpret_header(lines)
        
        # now interpret the ephemeris
        
        line = rinex_file.readline()
        while line != "" :
            lines = []
            lines.append(line)
            for ii in range(7) :
                line = rinex_file.readline()
                if line == "" :
                    break
                
                lines.append(line)
            
            if ii != 6 :
                break # end of file reached
            else :
                if lines[0][0] == 'G' :
                    eph = gps_rinex_ephemeris(lines)
                    
                    if eph not in eph_list :
                        eph_list.append(eph)
                    
                # Add here other ephemeris
            
                line = rinex_file.readline()
            
        # finally close the rinex file
        rinex_file.close()
        
        # return the ephemeris list
        return eph_list