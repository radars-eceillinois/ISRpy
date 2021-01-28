#! /usr/bin/env python
import numpy as np
import pylab as py
import datetime as dtm
import time as tmm
import os
import glob
from pdb import set_trace #useful for debugging
from scipy.constants import c
# History
# August 5, 2011 by Pablo M. Reyes
#   - modifications to handle filenames of the form *.pdata
#     done by using split() to get the filename without the path

shortH_t=np.dtype([
    ('length','<u4'),
    ('version','<u2'),
    ('datablock','<u4'),
    ('ltime','<u4'),
    ('milsec','<u2'),
    ('timezone','<i2'),
    ('dstflag','<i2'),
    ('errorCount','<u4')
    ])
sysH_t=np.dtype([
    ('nHeaderLength','<u4'),
    ('nNumSamples','<u4'),
    ('nNumProfiles','<u4'),
    ('nNumChannels','<u4'),
    ('nADCResolution','<u4'),
    ('nPCDIOBusWidth','<u4'),
    ])
rcfH_t=np.dtype([
    ('nHeaderLength','<u4'),
    ('nExpType','<u4'),
    ('nNTx','<u4'),
    ('fIpp','<f4'),
    ('fTxA','<f4'),
    ('fTxB','<f4'),
    ('nNumWindows','<u4'),
    ('nNumTaus','<u4'),
    ('nCodeType','<u4'),
    ('nLine6Function','<u4'),
    ('nLine5Function','<u4'),
    ('fClock','<f4'),
    ('nPrePulseBefore','<u4'),
    ('nPrePulseAfter','<u4'),
    ('sRangeIPP','<a20'),
    ('sRangeTxA','<a20'),
    ('sRangeTxB','<a20'),
    ])
procfH_t=np.dtype([
    ('nHeaderLength','<u4'),
    ('nDataType','<u4'),
    ('SizeOfDataBlock','<u4'),
    ('ProfilesperBlock','<u4'),
    ('DataBlocksperFile','<u4'),
    ('nNumWindows','<u4'),
    ('ProcessFlags','<u4'),
    ('CoherentIntegrations','<u4'),
    ('IncoherentIntegrations','<u4'),
    ('TotalSpectra','<u4')
    ])
samwin_t=np.dtype([
    ('h0','<f4'),
    ('dh','<f4'),
    ('nsa','<u4')])

class PROCFLAG:
    COHERENT_INTEGRATION=np.uint32(0x00000001)
    DECODE_DATA=np.uint32(0x00000002)
    SPECTRA_CALC=np.uint32(0x00000004)
    INCOHERENT_INTEGRATION=np.uint32(0x00000008)
    POST_COHERENT_INTEGRATION=np.uint32(0x00000010)
    SHIFT_FFT_DATA=np.uint32(0x00000020)
    DATATYPE_CHAR=np.uint32(0x00000040)
    DATATYPE_SHORT=np.uint32(0x00000080)
    DATATYPE_LONG=np.uint32(0x00000100)
    DATATYPE_INT64=np.uint32(0x00000200)
    DATATYPE_FLOAT=np.uint32(0x00000400)
    DATATYPE_DOUBLE=np.uint32(0x00000800)
    DATAARRANGE_CONTIGUOUS_CH=np.uint32(0x00001000)
    DATAARRANGE_CONTIGUOUS_H=np.uint32(0x00002000)
    DATAARRANGE_CONTIGUOUS_P=np.uint32(0x00004000)
    SAVE_CHANNELS_DC=np.uint32(0x00008000)
    DEFLIP_DATA=np.uint32(0x00010000)
    DEFINE_PROCESS_CODE=np.uint32(0x00020000)
    DATATYPE_MASK=np.uint32(0x00000FC0)

class PROCnDType:
    RAWDATA=np.uint32(0x00000000)
    SPECTRA=np.uint32(0x00000001)
    TWODIFF_PROCESS_WIN=np.uint32(0x00001000)
    THREEDIFF_PROCESS_WIN=np.uint32(0x00010000)
    FOURDIFF_PROCESS_WIN=np.uint32(0x00011000)
    SAVE_INCOH_INT_TIME_AVER=np.uint32(0x00100000)

class RCLnFun:
    NONE=0
    FLIP=1
    CODE=2
    SAMPLING=3
    LIN6DIV256=4
    SYNCHRO=5
class rcH:
    """will contain radar controller parameters of the JRO file header"""
    def __init__(self,fid):
        Flip1=0
        Flip2=0
        #read radar controller fixed header
        self.fixH=np.fromfile(fid,rcfH_t,1)
        #Start reading dynamic RC header
        #read RC Sampling Windows.
        self.Windows=np.fromfile(fid,samwin_t,self.fixH['nNumWindows'])
        self.Taus=np.fromfile(fid,'<f4',self.fixH['nNumTaus']) #read RC Taus.
        self.NumCodes=0
        self.NumBauds=0
        self.Codes=[]
        if self.fixH['nCodeType'] != 0:          #if a code was used
            self.NumCodes=np.fromfile(fid,'<u4',1)   #read number of codes
            self.NumBauds=np.fromfile(fid,'<u4',1)   #read number of Bauds
            self.Codes=np.empty([self.NumCodes,self.NumBauds],dtype='u1')
            for ic in range(self.NumCodes):
                #read 'ic' code
                tempc=np.fromfile(fid,'u1',4*np.ceil(self.NumBauds/32.))
                self.Codes[ic]=np.unpackbits(tempc[::-1])[-1*self.NumBauds:]
            #convert code from '0's and '1's to double +-'1'
            self.Codes=2.0*self.Codes-1.0
        if self.fixH['nLine5Function'] == RCLnFun.FLIP:
            self.rcFlip1=np.fromfile(fid,'<u4',1)  #read Flip1
        if self.fixH['nLine6Function'] == RCLnFun.FLIP:
            self.rcFlip2=np.fromfile(fid,'<u4',1)   #read Flip2

class procH:
    """will contain processing parameters of the JRO file header"""
    def __init__(self,fid):
        self.fixH=np.fromfile(fid,procfH_t,1)  #read processing fixed header
        #Start reading dynamic processing header
        #read proc Sampling Windows.
        self.Windows=np.fromfile(fid,samwin_t,self.fixH['nNumWindows'])
        #pdb.set_trace()
        #read spectra Combinations
        self.SpcCombinations = np.fromfile(fid, '<u1',
                2 * self.fixH['TotalSpectra'][0])
        self.SpcCombinations = self.SpcCombinations.reshape(
                self.fixH['TotalSpectra'][0], 2)
        self.NumCodes=0
        self.NumBauds=0
        self.Codes=[]
        if (self.fixH['ProcessFlags'] & PROCFLAG.DEFINE_PROCESS_CODE
                == PROCFLAG.DEFINE_PROCESS_CODE ):
            self.NumCodes=np.fromfile(fid,'<u4',1)   #read number of codes
            self.NumBauds=np.fromfile(fid,'<u4',1)   #read number of Bauds
            self.Codes=np.fromfile(fid, '<f4', NumCodes*NumBauds).reshape(
                    NumBauds, NumCodes)

class finfo:
    """will contain the JRO file headers & calculated experiment parameters"""
    def __init__(self,filename):
        fid=open(filename,'rb')
        self.fid=fid
        self.ThisFileSize=os.fstat(self.fid.fileno())[6]
        #read short fixed header into class
        self.shortH=np.fromfile(fid,shortH_t,1)
        self.sysH=np.fromfile(fid,sysH_t,1)	#read fixed system header
        self.rcH=rcH(fid)					#read RC header
        self.procH=procH(fid)				#read processing header
        self.fid.close()
        #correct header errors if needed
        self = correct_header_errors(self,filename)
        #add experiment parameters to class
        self = add_experiment_parameters(self)

    def getVals(self):
        """method to obtain all the values of an instantiated finfo"""
        #copy not to affect later with the del command
        pars=self.__dict__.copy()
        if pars.has_key('fid'):
            del pars['fid']
        return pars

def correct_header_errors(finfo,filename):
    #correcting known errors for JRO data
    fname = filename.split(os.sep).pop() #splits the path and pops up the fname
    yy = int(fname[1:5]) #year from filename
    dd = int(fname[5:5 + 3])  #doy from filename
    if (yy==2005 and dd in [74,75,76,105,115,116,117,164,165,166,167,168,248,
            249,250,251,346,347,348,354,355,356] or yy==2006 and dd in [93,94,
            95,96,213,214,215,248,249,250,338,339,340,341] or yy==2007 and
            dd in [170,171,172,173,174] ):
        #needed to correct wrong parameters in rawdata header
        if finfo.rcH.fixH['fIpp']==300.15 : finfo.rcH.fixH['fIpp']=300.0  #MST
        elif finfo.rcH.fixH['fIpp']==2997.0 : finfo.rcH.fixH['fIpp']=999.0 #ISR
        if finfo.rcH.fixH['fTxA']==45.0 : finfo.rcH.fixH['fTxA']=47.7
    return finfo

def add_experiment_parameters(finfo):
    """compute experiment parameters from JRO file headers and add to finfo"""
    finfo.lheader_length=finfo.shortH['length'][0]
    finfo.sheader_length=24
    setdigits = 3   #number of digits used for the set file set number
    fname = finfo.fid.name.split(os.sep).pop() #splits the path and gives fname
    finfo.year = int(fname[1:5])
    finfo.doy = int(fname[5:5+3])
    finfo.set = int(fname[5+3:5+3+setdigits])
    finfo.startime=finfo.shortH['ltime'][0]+finfo.shortH['milsec'][0]/1e3
    finfo.timezone_sec=finfo.shortH['timezone'][0]*60.
    #the information about the clock is in type float!!!
    finfo.dClock = (np.double(np.fix(finfo.rcH.fixH['fClock']*1000000.)[0])
                    / 1000000.)
    dResol=0.15/finfo.dClock #in kilometers
    finfo.ipp=dResol*np.round(finfo.rcH.fixH['fIpp'][0]/dResol)
    finfo.pw=dResol*np.round(finfo.rcH.fixH['fTxA'][0]/dResol)
    finfo.txa=dResol*np.round(finfo.rcH.fixH['fTxA'][0]/dResol)
    finfo.txb=dResol*np.round(finfo.rcH.fixH['fTxB'][0]/dResol)
    finfo.num_hei=sum(finfo.procH.Windows['nsa'])
    finfo.num_win=finfo.procH.fixH['nNumWindows'][0]
    finfo.first_heigth=dResol*np.round(finfo.procH.Windows['h0'][0]/dResol)
    finfo.spacing=dResol*np.round(finfo.procH.Windows['dh'][0]/dResol)
    finfo.samples_win=finfo.procH.Windows['nsa'][0]
    if finfo.procH.fixH['nDataType'][0] & PROCnDType.SPECTRA == \
            PROCnDType.SPECTRA:
        finfo.num_chan=sum(finfo.procH.SpcCombinations[:,0] ==
            finfo.procH.SpcCombinations[:,1])
        finfo.num_pairs=finfo.procH.fixH['TotalSpectra'][0] - finfo.num_chan
    else :
        finfo.num_chan=finfo.sysH['nNumChannels'][0]
        finfo.num_pairs=0
    finfo.num_prof=finfo.procH.fixH['ProfilesperBlock'][0]
    finfo.num_coh=max(finfo.procH.fixH['CoherentIntegrations'][0],1)
    finfo.num_incoh=max(finfo.procH.fixH['IncoherentIntegrations'][0],1)
    finfo.bytes_block=(finfo.procH.fixH['SizeOfDataBlock'][0] +
                        finfo.sheader_length)
    finfo.blocks_file=finfo.procH.fixH['DataBlocksperFile'][0]
    finfo.ThisFile_blks=(np.round((finfo.ThisFileSize-finfo.lheader_length +
                        finfo.sheader_length)/finfo.bytes_block))
    dtmask = finfo.procH.fixH['ProcessFlags'][0] & PROCFLAG.DATATYPE_MASK
    if dtmask == 0 :
        finfo.data_type = 4
        print "error reading process flags, data type used instead: float"
    else :
        finfo.data_type=int(np.log2(dtmask) - np.log2(PROCFLAG.DATATYPE_CHAR))
    finfo.taus=finfo.rcH.Taus
    finfo.hrange=np.arange(finfo.num_hei)*finfo.spacing+finfo.first_heigth
    #useful for decoding
    finfo.code=finfo.rcH.Codes
    #PROC information
    finfo.ProcInfo=ProcInfo(finfo.procH.fixH)
    #extra useful parameters for decoding
    if finfo.code!=[]:
        finfo.subbauds=np.round(finfo.txa / finfo.spacing / np.double(
                finfo.rcH.NumBauds))
        finfo.subcode=expandcode(finfo.rcH.Codes,finfo.subbauds)
        finfo.num_codes=finfo.subcode.shape[0]
        finfo.num_bauds=finfo.subcode.shape[1]
        if finfo.ProcInfo.decoded:
            finfo.deco_num_hei=finfo.num_hei
            finfo.deco_hrange=finfo.hrange
        else:
            finfo.deco_num_hei=finfo.num_hei-finfo.num_bauds+1
            finfo.deco_hrange=finfo.hrange[:-finfo.num_bauds+1]
    #extra useful parameters to obtain velocities
    #this frequency is used to have a lambda closer to 6m
    finfo.rad_freq=49.92e6
    finfo.rad_lambda=c/finfo.rad_freq
    finfo.ipp_secs=2.*finfo.ipp*1000.*finfo.num_coh/c
    finfo.freq_interval=1./finfo.ipp_secs;
    finfo.velrange=finfo.freq_interval*finfo.rad_lambda/2.
#            profs=self.Hinfo.num_prof
#            self.vels=(np.arange(profs)-profs/2.)/profs*self.velrange
#    finfo.fft_code=np.fft.fft(Hinfo.subcode,Hinfo.num_hei,1).conj()

    return finfo

class ProcInfo:
    """will contain information of the processes done at acquisition time"""
    def __init__(self,procHfixH):
        if (procHfixH['ProcessFlags'] & PROCFLAG.DECODE_DATA
                == PROCFLAG.DECODE_DATA):
            self.decoded=True
        else: self.decoded=False
        if (procHfixH['ProcessFlags'] & PROCFLAG.SPECTRA_CALC
                == PROCFLAG.SPECTRA_CALC):
            self.spectra=True
        else: self.spectra=False
        if (procHfixH['ProcessFlags'] & PROCFLAG.INCOHERENT_INTEGRATION
                == PROCFLAG.INCOHERENT_INTEGRATION):
            self.incoh_int=True
        else: self.incoh_int=False
        if (procHfixH['ProcessFlags'] & PROCFLAG.DEFLIP_DATA
                == PROCFLAG.DEFLIP_DATA):
            self.deflip=True
        else: self.deflip=False
        if (procHfixH['ProcessFlags'] & PROCFLAG.DATAARRANGE_CONTIGUOUS_CH
                == PROCFLAG.DATAARRANGE_CONTIGUOUS_CH):
            self.data_contiguous_ch=True
        else: self.data_contiguous_ch=False

def expandcode(code,expnum):
    """tx codes are expanded by the # of samples per baud to produce subcode"""
    return np.reshape(np.tile(np.reshape(code,(np.size(code),1)),(1,expnum)),
        (code.shape[0],code.shape[1]*expnum))

def readH_continue(fid):
    """move past the next header and return the short header"""
    try:
        shortH=np.fromfile(fid,shortH_t,1)
        if shortH.size>0:
            if shortH['length']>24 :
                fid.seek(shortH['length']-24,1)
            elif shortH['length']<24:
                return []
        else:
            return []
        return shortH
    except:
        return []

def readRawBlk(finfo,profs=None):
    """return block of complex raw data of specified # of profiles or all"""
     #0:Int8, 1:Int16, 2:Int32, 3:Int64, 4:Float, 5: Double
    dtyps=['<i1','<i2','<i4','<i8','<f4','<f8']
    dtyp=dtyps[finfo.data_type]
    dt=np.dtype([('real',dtyp),('imag',dtyp)])

    if profs==None: profs=finfo.num_prof
    pts2read=profs*finfo.num_chan*finfo.num_hei
    try:
        temp=np.fromfile(finfo.fid,dt,pts2read)
    except:
        return []
    if temp.nbytes==0:
        return []
    temp=temp.reshape((profs,finfo.num_hei,finfo.num_chan))
    return temp['imag']*1j+temp['real']

def readProcBlk(finfo,profs=None):
    """return block of proccess data of specified # of profiles or all"""
    #0:Int8, 1:Int16, 2:Int32, 3:Int64, 4:Float, 5: Double
    dtyps=['<i1','<i2','<i4','<i8','<f4','<f8']
    dtyp=dtyps[finfo.data_type]
    dtcmplx=np.dtype([('real',dtyp),('imag',dtyp)])
    dtreal=np.dtype(dtyp)
    if profs==None: profs=finfo.num_prof
    #reading the self spectra
    pts2read=profs*finfo.num_chan*finfo.num_hei
    try:
        tempspc=np.fromfile(finfo.fid,dtreal,pts2read)
    except:
        return [],[],[]
    if finfo.num_pairs > 0: #if cross spectra
        #reading the cross spectra
        pts2read=profs*finfo.num_pairs*finfo.num_hei
        try:
            tempcspc=np.fromfile(finfo.fid,dtcmplx,pts2read)
        except:
            return [],[],[]
    if finfo.procH.fixH['ProcessFlags'][0] & PROCFLAG.SAVE_CHANNELS_DC == \
                PROCFLAG.SAVE_CHANNELS_DC:
        #reading the DC
        pts2read=finfo.num_chan*finfo.num_hei
        try:
            tempDC=np.fromfile(finfo.fid,dtcmplx,pts2read)
        except:
            return [],[],[]

    if tempspc.nbytes==0 or tempcspc.nbytes==0 or tempDC.nbytes==0:
        return [],[],[]

    tempspc=tempspc.reshape((finfo.num_chan,finfo.num_hei,profs))
    tempcspc=tempcspc.reshape((finfo.num_pairs,finfo.num_hei,profs))
    tempDC=tempDC.reshape((finfo.num_chan,finfo.num_hei))
    fcspc = tempcspc['imag']*1j+tempcspc['real']
    fDC = tempDC['imag']*1j+tempDC['real']
    return (tempspc, fcspc, fDC)

def skip_blocks(finfo,blks2skip):
    """skip blks2skip header+data blocks from the beginning of the datafile"""
    dbyts=[1,2,4,8,4,8]
    databytes=dbyts[finfo.data_type]
    if finfo.procH.fixH['nDataType'][0] & PROCnDType.SPECTRA == \
            PROCnDType.SPECTRA:
        #self spectra
        bytesperblk=databytes*finfo.num_prof*finfo.num_chan*finfo.num_hei
        if finfo.num_pairs > 0: # if cross spectra calculated
            #cross spectra (real and imaginary)
            bytesperblk = (bytesperblk +
                2 * databytes * finfo.num_prof * finfo.num_pairs*finfo.num_hei)
        if finfo.procH.fixH['ProcessFlags'][0] & PROCFLAG.SAVE_CHANNELS_DC == \
                PROCFLAG.SAVE_CHANNELS_DC:
            # complex DC value for each channel, and height
            bytesperblk = (bytesperblk +
                2 * databytes * finfo.num_chan*finfo.num_hei)
    else:
        #rawdata (real and imaginary)
        bytesperblk = 2* databytes*finfo.num_prof*finfo.num_chan*finfo.num_hei
    print bytesperblk
    print finfo.procH.fixH['SizeOfDataBlock'][0]
    if finfo.procH.fixH['SizeOfDataBlock'][0] != bytesperblk:
        print 'error calculating size of the block'
        bytesperblk = finfo.procH.fixH['SizeOfDataBlock'][0]
    bytes2skip=(finfo.lheader_length-24+blks2skip*(bytesperblk+24))
    try:
        #0: from the beginning of the file,1: current position,2:end
        finfo.fid.seek(bytes2skip,0)
    except:
        return -1   # returns -1 if there was an i/o error
    if finfo.fid.tell()!=bytes2skip:
        return -2   #returns -2 if the skip bytes does not coincide
    else:
        #returns the number of bytes skiped if successfully done
        return bytes2skip

def create_finf_from_savemat(matobj):
    class finf: #initiating class finf in order to populate it later
        """
        This class will contain members similar to those of a finfo class
        extracting values created by the function getVals(self) and saved
        as a matlab file. """
        pass
    class emptyc:
        "This is a generic class for rcH and procH"
        pass
    for keynames0 in matobj.dtype.names:
        obj = matobj[keynames0]
        while obj.dtype == 'object': #eliminate the 'object' type
            obj = obj[0]
        if obj.dtype.names == None : #None means single value or array
            if obj.size == 1: #for single 1 element values
                finf.__dict__[keynames0] = obj[0][0]
            else: #for empty values or arrays
                finf.__dict__[keynames0] = obj.squeeze()
        elif obj.dtype.names[1]=='nNumSamples':#sysH_t dtype
            #different thatn None means a structure of names
            finf.__dict__[keynames0] = obj.astype(sysH_t).squeeze()
        elif obj.dtype.names[1]=='version':#shortH_t dtype
            finf.__dict__[keynames0] = obj.astype(shortH_t).squeeze()
        else:
            class1 = finf.__dict__[keynames0] = emptyc() #starting class
            for keyname1 in obj.dtype.names:
                obj1 = obj[keyname1]
                while obj1.dtype == 'object': #eliminate the 'object' type
                    obj1 = obj1[0]
                if obj1.dtype.names == None:
                    if obj1.size == 1: #for single 1 element values
                        class1.__dict__[keyname1] = obj1[0][0]
                    else: #for empty values or arrays
                        class1.__dict__[keyname1] = obj1.squeeze()
                elif obj1.dtype.names[1]=='dh':#samwin_t dtype
                    # a structure inside another structure
                    class1.__dict__[keyname1] = obj1.astype(samwin_t).squeeze()
                elif obj1.dtype.names[1]=='nDataType':#procfH_t dtype
                    class1.__dict__[keyname1] = obj1.astype(procfH_t).squeeze()
                elif obj1.dtype.names[1]=='nExpType':#rcfH_t dtype
                    temp = class1.__dict__[keyname1] = np.array(
                            '', dtype = rcfH_t)
                    for tname in rcfH_t.names:
                        if tname[:6] == 'sRange':
                            temp[tname] = obj1[tname][0][0].tostring(
                                    ).replace('\x00','')
                        else:
                            temp[tname] = obj1[tname][0][0].squeeze()
    return finf

def runTest():
    """
    import jropack.jroread as jro
    [finf,data]=jro.runTest()
    """
    #import pylab as py
    fig = py.figure()	# initiate graphic backend
    py.rcParams.update({'xtick.labelsize':10,'ytick.labelsize':10})
    foldername = 'd2009013'
    #foldername = '../../data/MSTISR/rawdata/y2009/d2009018/ISR'
    fullname=foldername+'/D2009013007.r'
    #fullname=foldername+'/D2009018125.r'
    finf=finfo(fullname) #read header and close file ready for further action
    finf.fid=open(fullname,'rb')		#reopen file
    shortH=readH_continue(finf.fid)	#move past first header
    data=readRawBlk(finf)			#read data block
    #py.plot(data.real[0,:,:],data.imag[0,:,:],'or')
    py.plot(finf.hrange,data.real[0,:,:])
    py.xlabel('Height (km)')
    #shortH=readH_continue(finf.fid)	#move past next header
    fig.show()
    finf.fid.close()
    py.show()
    return [finf,data]

if __name__=='__main__':
	[finf,data]=runTest()
