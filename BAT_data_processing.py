from astropy.io import fits
import os
import numpy as np
from math import floor
from subprocess import run


class BAT_data(object):
    """
    This object makes BAT fits-files and applies correction for the spacecraft slew to the DRM
    """
    
    def __init__(self, GRB_name=None, trig_ID=None, dT0=0, MET=None, ToF=0, slew_interval=None, sp_interval=None, \
                 sp_ID=None, path2BATdata=None, new_folder=None, verbose=False):
        
        """
        Initialisation of the BAT_data object.
        
        Parameter:
            - GRB_name: Burst GCN name
            - trigID: Long BAT trigger ID
            - dT0: Difference in KW and BAT trigger times (s)
            - MET: Swift mission elapsed time (s)
            - ToF: 'Time of flight' - Time of light propagation between KW and Swift (s)
            - slew_interval: Time interval when the spacecraft slews (s)
            - sp_interval: Spectrum time interval (s; relative KW trigger)
            - spID: Designation of the spectrum
            - path2BATdata: Path to the directory containing the BAT folder 'trigID-results/'
            - new_folder: Folder where the new BAT fits-files will be collected
        """
        
        self.verbose = verbose
    
        if not GRB_name:
            raise TypeError('Burst name should be given')
        elif not isinstance(GRB_name, str):
            raise TypeError('Burst name should be a string')
        else:
            self._name = GRB_name

        if not isinstance(sp_ID, str):
            raise TypeError('Spectrum designation should be given')
        else:
            self._idsp = sp_ID    

        self._dT0 = float(dT0)
        self._ToF = float(ToF)
    
        if MET is None:
            raise TypeError('MET should be given')
        else:
            self._MET = float(MET)

        if slew_interval is None:
            raise TypeError('Slew time interval should be given')
        else:
            if not isinstance(slew_interval, list):
                raise TypeError('Slew time interval should be a list')
            else:
                self._slew_interval = list(map(float, map('{:.3f}'.format, slew_interval)))

        if sp_interval is None:
            raise TypeError('Spectrum time interval should be given')
        else:
            if not isinstance(sp_interval, list):
                raise TypeError('Spectrum time interval should be a list')
            else:
                self._sp_interval = sp_interval
            
        if not trig_ID:
            raise TypeError('Trigger ID should be given')
        elif not isinstance(trig_ID, str):
            raise TypeError('Trigger ID should be a string')
        elif len(trig_ID) != 11:
            raise TypeError('Trigger ID should be an eleven-digit string')
        else:
            self._trig = trig_ID

        if not path2BATdata:
            raise TypeError('Path to the directory containing the BAT folder trigID-results/ should be given')
        else:
            self._path2BAT = path2BATdata

        if not new_folder:
            raise TypeError('Folder where the new BAT fits-files will be collected should be given')
        else:
            self._new_folder = new_folder

        if self.verbose:
            print('Burst name:', self._name)
            print('Spectrum:', self._idsp, 'lasts from', self._sp_interval[0], 's to', self._sp_interval[1], 's')
            print('Long BAT trigger ID:', self._trig)
            print('Swift Mission Elapsed Time:', self._MET, 's')
            print('Time of light propagation between KW and Swift:', self._ToF, 's')
            print('Difference in KW and BAT trigger times:', self._dT0, 's')
            print('Time interval when the spacecraft slews: from', self._slew_interval[0], 's to', self._slew_interval[1], 's')
            print('Path to the BAT data folder:', self._path2BAT)
            print('Folder where the new BAT fits-files will be collected:', self._new_folder)
        

    def make_rsp(self):
        """
        Calculating DRM weights and running a slew-corected DRM
        """
                
        if self.verbose:
            print(self.make_rsp.__doc__)
        
        slice_times = list()
        slice_cnts = list()
        cnts = dict()
        
        t_i, t_f = self._sp_interval
        inslew_i, inslew_f = self._slew_interval
        
        if t_i <= inslew_i and inslew_f <= t_f:
            t_slew = inslew_f - inslew_i

            if t_slew >= 6:
                n = floor(t_slew / 5)
                if t_slew % 5 >= 1: n += 1

                slice_cnts.append(get_bat_counts(t_i, inslew_i))
                slice_times.append([t_i, inslew_i])

                for i in range(0, n):
                    t1 = t_i + i*5
                    t2 = t_i + 1 + i*5 if t2 <= inslew_f else inslew_f
                    slice_cnts.append(get_bat_counts(t1, t2))
                    slice_times.append([t1, t2])

                slice_cnts.append(get_bat_counts(inslew_f, t_f))
                slice_times.append([inslew_f, t_f])

                if self.verbose:
                    print(f'GRB {self._name}: DRM will be corrected for the first type s/c slew')

        elif t_i <= inslew_i and t_f <= inslew_f:
            t_slew = t_f - inslew_i

            if t_slew >= 6 :
                    
                n = floor(t_slew / 5)
                if t_slew % 5 >= 1: n += 1  

                slice_cnts.append(get_bat_counts(t_i, inslew_i))
                slice_times.append([t_i, inslew_i])

                for i in range(0, n):
                        
                    t1 = inslew_i + i*5
                    t2 = float()
                    t2 = inslew_i + 1 + i*5 if t2 <= t_f else t_f
   
                    slice_cnts.append(get_bat_counts(t1, t2))
                    slice_times.append([t1, t2])

                if self.verbose:
                    print(f'GRB {self._name}: DRM will be corrected for the second type s/c slew')

        elif inslew_i <= t_i and inslew_f <= t_f:
            t_slew = inslew_f - t_i

            if t_slew >= 6:
                n = floor(t_slew / 5)
                if t_slew % 5 >= 1: n += 1

                for i in range(0, n):
                        
                    t1 = t_i + i*5
                    t2 = t_i + 1 + i*5 if t2 <= inslew_f else inslew_f

                    slice_cnts.append(get_bat_counts(t1, t2))
                    slice_times.append([t1, t2])

                slice_cnts.append(get_bat_counts(inslew_f, t_f))
                slice_times.append([inslew_f, t_f])

            if self.verbose:
                print(f'GRB {self._name}: DRM will be corrected for the third type s/c slew')

        elif inslew_i <= t_i and t_f <= inslew_f:
            t_slew = t_f - t_i

            if t_slew >= 6:
                n = floor(t_slew / 5)
                if t_slew % 5 >= 1: n += 1

                for i in range(0, n):
                        
                    t1 = t_i + i*5
                    t2 = t_i + 1 + i*5 if t2 <= t_f else t_f

                    slice_cnts.append(get_bat_counts(t1, t2))
                    slice_times.append([t1, t2])

            if self.verbose:
                print(f'GRB {self._name}: DRM will be corrected for the fourth type s/c slew')

        total_cnts = np.sum(slice_cnts)
        weights = slice_cnts / total_cnts

        if len(slice_times) == 0:
            slice_times = [[t_i, t_f]]
            weights = [1.]
            
        if self.verbose:
            print('Time intervals for aux DRMs:', *slice_times)
            print('Weights of the aux DRMs:', *weights)

        self.make_new_bat_drm(slice_times, weights)

            
    def get_bat_counts(self, t1, t2):
        """
        Extracting counts for a certain time interval from a BAT light curve
        """
        
        if self.verbose:
            print(self.get_bat_counts.__doc__)
            
        if t1 >= t2:
            raise TypeError('t2 should be greater than t1')

        fits_filename = self._path2BAT + f"/{self._trigID}-results/lc/sw{self._trigID}b_4chan_64ms.lc"
    
        with fits.open(fits_filename) as hdul:
            data = hdul[1].data
            data = np.array([*data])

            MET = hdul[1].header['TRIGTIME']
            data[:,0] = data[:,0] - MET - dT0 - ToF

            cnts = data[t1 < data[:,0]]
            cnts = cnts[cnts[:,0] < t2]

            sum_cnts = np.sum(np.sum(cnts[:,1]))
            return 0.064 * sum_cnts
      
        
    def make_new_bat_drm(self, slice_times, weights):
        """
        Making a slew-corrected DRM
        """
        
        if self.verbose:
            print(self.make_new_bat_drm.__doc__)

        evt_file = self._path2BAT + f"/{self._name}/{self._trig}-results/events/sw{self._trig}b_all.evt"
        mask_file= self._path2BAT + f"/{self._name}/{self._trig}-results/auxil/sw{self._trig}b_qmap.fits"
        ray_tracer_file = self._path2BAT + f"/{self._name}/{self._trig}-results/auxil/sw{self._trig}b_all.evaux"
        rmf_file = self._path2BAT + f"/{self._name}/{self._trig}-results/{self._new_folder}/BAT_{self._trig}_sp{self._idsp}_averaged.rmf"
        rmf_list = list()
        aux_list = list()

        for times in slice_times:

            pha_tstart, pha_tstop = times + (self._MET + self._ToF + self._dT0)* np.ones(2)

            gtinum=1

            pha_file = self._path2BAT + f"/{self._name}/{self._trig}-results/{self._new_folder}/sw{self._trig}b_{gtinum}_aux.pha"
        
            # Extract spectra rebinned in time 
            sys = f"batbinevt infile={evt_file} outfile={pha_file} outtype=PHA timedel=0 timebinalg=u \
            energybins=CALDB:80 tstart={pha_tstart} tstop={pha_tstop} outunits=RATE detmask={mask_file} \
            clobber=yes ecol=ENERGY weighted=YES"
            if self.verbose:
                print(sys)
            run(sys, shell=True)
            
            aux_list.append(pha_file)
        
            # Update BAT ray tracing columns in spectral files 
            sys = f"batupdatephakw infile={pha_file} auxfile={ray_tracer_file}"
            if self.verbose:
                print(sys)
            run(sys, shell=True)

            rsp_file = self._path2BAT + f"/{self._name}/{self._trig}-results/{self._new_folder}/sw{self._trig}b_{gtinum}_aux.rsp"

            sys = f"batdrmgen infile={pha_file} outfile={rsp_file} hkfile=NONE clobber=yes"
            if self.verbose:
                print(sys)
            run(sys, shell=True)

            rmf_list.append(rsp_file)
            rmf_line = ','.join(rmf_list)
            weight_line = ','.join(map(str, weights))

            gtinum += 1 
            

            sys = f"addrmf {rmf_line} {weight_line} rmffile={rmf_file}"
            if self.verbose:
                print(sys)
            run(sys, shell=True)

            for rmf in rmf_list:
                os.remove(rmf)

            for aux in aux_list:
                os.remove(aux)
                
    def make_pha(self):
        """
        Making a pha file
        """
        
        if self.verbose:
            print(self.make_pha.__doc__)

        pha_tstart, pha_tstop = self._sp_interval + (self._MET + self._ToF + self._dT0)* np.ones(2)
        
        evt_file = self._path2BAT + f"/{self._name}/{self._trig}-results/events/sw{self._trig}b_all.evt"
        pha_file = self._path2BAT + f"/{self._name}/{self._trig}-results/{self._new_folder}/BAT_{self._trig}_sp{self._idsp}.pha"
        mask_file = self._path2BAT + f"/{self._name}/{self._trig}-results/auxil/sw{self._trig}b_qmap.fits"
        ray_tracer_file = self._path2BAT + f"/{self._name}/{self._trig}-results/auxil/sw{self._trig}b_all.evaux"
        
        # Accumulate BAT event and DPH data into spectra, lightcurves or images:
        # Extract spectra rebinned in time 
        
        sys = f"batbinevt infile={evt_file} outfile={pha_file} outtype=PHA timedel=0.0 timebinalg=u \
        tstart={pha_tstart} tstop={pha_tstop} energybins=CALDB:80 outunits=RATE detmask={mask_file} clobber=YES \
        ecol=ENERGY weighted=YES"
        if self.verbose:
            print(sys)
        run(sys, shell=True)    
        
        # Update BAT ray tracing columns in spectral files 
        sys = f"batupdatephakw infile={pha_file} auxfile={ray_tracer_file}"
        if self.verbose:
            print(sys)
        run(sys, shell=True)
        
        # Apply BAT spectral systematic error vector
        sys = f"batphasyserr infile={pha_file} syserrfile=CALDB"
        if self.verbose:
            print(sys)
        run(sys, shell=True)
