from numpy import *
from struct import *
import datetime

class JROfileread: # create a class instead of a function if you need to rememeber info in between calls

    """
    instantiate/open a JRO data file by its name using the class name and
    read the file block by block using the class method nextdatablock:

    Example:
    foo=JROfileread(filename)
    data=foo.nextdatablock() --- which is a 3D complex voltage array of [cycles,samples,channels] dimensions

    also explore class parameters and methods by tabbing "foo."" in ipython or jupyter
    """

    def __init__(self,filename): # instantiate by filename input ... self variable always needed in a class
        fid=open(filename,'rb')
        self.fid=fid # in a class self.xxx variables are generated for remembering between calls

    def nextdatablock(self): # method to read headaer and data block by block
        fid=self.fid # remembering something memorized
        shortH=fid.read(24) # read short header bytes, unpack, and go to the data block
        hlength,version,datablock,ltime,msec,timezone,dstflag,errorCount=unpack('<ihiihhhi',shortH)

        if hlength>24: # but if first block, need to read other parameters from extra headers before the data block
            sysH=fid.read(24) # system header
            length,nSamples,nPulses,nChannels,nADCResolution,nPCDIOBusWidth=unpack('<iiiiii',sysH)
            self.nSamples=nSamples # memorizing useful stuff, for inside and outside use later on!!!
            self.nChannels=nChannels
            # stored self.xxx variables can be accessed as foo.xxx from outside

            rcH=fid.read(56) # radar controller header
            length,Exp,nCycles,IPP,txA,txB,nWin,nTaus,CodeType,n6,n5,clock,ppb,ppa=unpack('<iiifffiiiiifii',rcH)
            self.nCycles=nCycles # more useful stuff
            self.IPP=IPP

            fid.seek(60,1) # skip some useless stuff

            sampleH=fid.read(12) # sampling header
            h0,dh,nsa=unpack('<ffi',sampleH)
            self.h0=h0
            self.dh=dh
            self.hts=h0+dh*arange(nsa)

            fid.seek(hlength,0) # now go to data block

        self.datablock=datablock # read data block after here
        self.starttime=ltime+msec/1000.
        self.date=datetime.datetime.fromtimestamp(ltime).strftime('%Y-%m-%d %H:%M:%S')

        allsamples=self.nSamples*self.nChannels*2
        data2read=self.nCycles*allsamples
        data=fid.read(data2read*2)  # read and unpack data block
        data=array(unpack('h'*data2read,data))

        data=reshape(data,(self.nCycles,self.nSamples,self.nChannels,2))
        data=data[:,:,:,0]+1j*data[:,:,:,1] # complex voltage array of [cycles,samples,channels] dimensions

        return data
