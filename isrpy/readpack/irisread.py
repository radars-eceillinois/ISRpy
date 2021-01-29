from os import path,stat
from scipy import constants
from numpy import dtype,fromfile,median,arange,array,argmin,complex64,ones,nan
from matplotlib.mlab import find
from time import gmtime
from calendar import timegm
class fileinfo:
    """    fileinfo initialization:

    from IRISread import fileinfo
    rawinfo = fileinfo("some filepath")

    Examples reading rawdata pulses
    1) Read the first 200 pulses saved in a rawdata file:
    raw = rawinfo.readNpulses(200)
    2) Read 200 pulses after skipping the first 100 pulses:
    raw = rawinfo.readNpulses(200,100)
    3) Read 200 pulses after skipping the first 100 pulses and return number of read pulses:
    raw,pulsesread = rawinfo.readNpulses(200,100,out_pulsesread=True)

    """
    def __init__(self,filepath):
        self.filepath = filepath
        timename = path.splitext(path.split(filepath)[1])[0]
        year,mon,day,HH,MM,SS = array(timename.split('.'),dtype=int)
        self.localtime_seconds = timegm((year,mon,day,HH,MM,SS,0,0,0))
        self.localtime_struct = gmtime(self.localtime_seconds)
        self.dtyp = dtype([('real','<f4'),('imag','<f4')])
        self.samp_rate = 400e3                  # Rx Sampling rate (Hz)
        self.dh = 1/self.samp_rate             #=2.5e-6 Height res. in seconds
        self.dhkm = 3e8/2/self.samp_rate/1000. # = 0.375 km. Approx height res
        # Exact height resolution in km :
        self.dhkm_exact = constants.c/2/self.samp_rate/1000. # = 0.3747405725 km
        self.dbytes = 4. # number of bytes
        self.IQ = 2.                 # In-Phase - In-Quadrature
        self.RXs = 2.                # number of receivers (Antennas)

        trysamp2read = 10 * self.RXs * 3200 # at least 10 pulses of largest case
        if not path.exists(filepath):
            print "File not found:",filepath
            return
        try:
            fid = open(filepath,'rb')
            pwr = fromfile(fid,self.dtyp,int(trysamp2read))
            pwr = pwr['real'] **2 + pwr['imag'] ** 2
            pwr = pwr.reshape(trysamp2read/self.RXs,self.RXs)
            fid.close()
        except:
            print "Error reading file:",filepath
            return
        self.nhts,self.h_offset = self.__GetIPPhts__(pwr[:,0])
        del(pwr)
        #ch0:less signal power, therefore is more desirable
        self.hts = arange(self.nhts) * self.dhkm
        self.IPP = self.dh * self.nhts # 8ms, 4ms, 5ms
        self.prf = 1./self.IPP   # pulse repetition frequency: 125 Hz
        self.IPPkm = self.IPP * 150e3 # 1200 km
        self.IPPkm_exact = self.IPP*constants.c/2000 # 1199.169832 km
        self.saved_pulses = int(stat(filepath).st_size / self.dbytes / self.IQ
                / self.RXs / self.nhts)

        self.frequency_span = 1. / self.IPP # Total span in Hz due to samp_rate
        self.radar_oper_freq = 49.8e6 # Hz
        self.radar_oper_lambda = constants.c /  self.radar_oper_freq # m
        self.Bragg_lambda = self.radar_oper_lambda / 2 #m
        self.vel_span = self.Bragg_lambda * self.frequency_span

    def vel(self,npulses):
        return (arange(npulses)-npulses/2.) / npulses * self.vel_span

    def __GetIPPhts__(self,pwr):
        peaklocations = find(pwr > pwr.max()*0.9)
        differences = peaklocations[1:] - peaklocations[:-1]
        median_diff = median(differences)
        return int(median_diff),peaklocations[0]

    def printParameters(self):
        for param in vars(self).keys():
            units = ''
            if param in ['IPPkm','IPPkm_exact','dhkm','dhkm_exact','hts']:
                units = 'km'
            elif param in ['IPP','dh']:
                units = 's'
            elif param in ['prf','frequency_span','samp_rate','radar_oper_freq']:
                units = 'Hz'
            elif param in ['Bragg_lambda','radar_oper_lambda']:
                units = 'm'
            if param == 'hts':
                print param,'= [%.3f, %.3f, %.3f, ..., %.3f, %.3f, %.3f]'%tuple(
                        self.hts[[0,1,2,-3,-2,-1]]),units
            else:
                print param,'=', getattr(self,param),units


    def readNpulses(self,pulses2read,pulses2skip=0,out_pulsesread=False):
        samp2skip = pulses2skip * self.RXs * self.nhts
        samp2read = pulses2read * self.RXs * self.nhts
        fid = open(self.filepath,'rb')
        fid.seek(samp2skip * self.IQ * self.dbytes,0)
        raw = fromfile(fid,self.dtyp,int(samp2read))
        pulsesread = raw.size / self.RXs / self.nhts
        fid.close()
        raw = raw['real'] + 1j * raw['imag']
        raw = raw.reshape(pulses2read,self.nhts,self.RXs)
        if out_pulsesread:
            return raw,pulsesread
        else:
            return raw

    def read_nIPP_Hkm(self,nIPP,Hkm,nhts=1):
        hindex = argmin(abs(Hkm-self.hts))
        fid = open(self.filepath,'rb')
        fid.seek((self.h_offset+hindex) * self.RXs * self.IQ * self.dbytes,1)
        raw = complex64(nan) * ones((nIPP,nhts,self.RXs),dtype=complex64)
        for i in range(int(nIPP)):
            raw[i,:,:] = fromfile(fid, dtype=(complex64,(self.RXs)), count=nhts )
            fid.seek(self.IQ*self.dbytes*self.RXs*(self.nhts-nhts),1)
        fid.close()
        return raw.squeeze()


