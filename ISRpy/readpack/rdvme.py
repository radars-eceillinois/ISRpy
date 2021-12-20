"""
Classes for reading records from the vme crate.


class rdvme---Reading aeronomy records from the VME system.
A class is defined; creating an instance requires an argument which is a string
containing the name of the file.
Assume the instance name is df, then operations such as seek and tell are 
available as df.f.tell() and df.f.seek(..., ...)
Methods:
 1. read a record,
 2. advance n records,
 3. go back to the beginning of the file
 4. move the file pointer to the next (often the first) occurrence of 'hdr_'
    (the mark for the beginning of a header)
Methods, where appropriate, return 1 for normal completion and 0 for EOF.
The header is read in a general way, using the info in the standard header to
identify the locations of the sps, ri, and pr subheaders.
The headers are decoded, and there is a tuple for each as data in the class instance
(that is, using df as above, df.stdhdr, df.spshdr, df.rihdr, and df.prhdr.)
Commonly used parameters also are in df, as int, float, or string, as appropriate.
Use dir(df), or however the instance is named, to get a list of instance contents.

class rdvmew---Like rdvme, except when EOF is found, this version of rdvmerec waits, and then checks each second to see if more has been written.
This class inherits from class rdvme, overriding __init__ and rdvmerec.  rdvmerec
returns 1, always, for compatibility.

class rdvmewf---Like rdvmew, but uses a file base and number. It opens and reads from
the next file in the numbered sequence of files when it is available.  Until it is
available, it continues try to read from the current file.
This class inherits from rdvmew, overriding __init__ and rdvmerec.  rdvmerec
returns 1, always, for compatibility.

No matter which class is used, the routine that reads the record is called
rdvemrec.  Just use the correct class for the job.

rdvmew and rdvmewf have test routines at the end of this file.
"""
from struct import *
import numpy as np
import time


# Class for reading records written by the VME data taking system
class rdvme:
    def __init__(self, filestring):
        self.f = open(filestring, 'rb')
        self.filestring = filestring

# method for reading a vme data taking record
# It knows about the standard data taking programs.
# (pwr, mracf, clp, topsd, dspc; these have floating point data output;
#  and rawdat; it has int16 output that is converted to float32)
# Float is read if the program hdr exists (length > 0)
# Otherwise int16 is read
# Returns 0 when EOF is encountered; normal return value is 1
    def rdvmerec(self):
        loc = self.f.tell()
        self.infostr = self.f.read(12)
        if len(self.infostr) == 0: return 0
        self.hdrinfo = unpack_from('>4sii', self.infostr)
        self.f.seek(loc, 0)
        self.hstr = self.f.read(self.hdrinfo[1])
        if len(self.hstr) == 0: return 0
        self.stdstr = '>4s2i8s7i16B16i'
        self.stdhdr = unpack_from(self.stdstr, self.hstr)
        self.ptype = self.stdhdr[3]
        self.ristr =  '>' + str(self.stdhdr[12]*16) + 'x' + str(self.stdhdr[14]//4) + 'i'
        self.prstr =  '>' + str(self.stdhdr[20]*16) + 'x' + str(self.stdhdr[22]//4) + 'i'
        self.spsstr = '>' + str(self.stdhdr[24]*16 + self.stdhdr[25]//16) + 'x8s7fi21f3sxi3f2if2if2i'
        self.rihdr =  unpack_from(self.ristr, self.hstr)
        if self.stdhdr[22] > 0: 
            self.prhdr = unpack_from(self.prstr, self.hstr)
        self.spshdr =  unpack_from(self.spsstr, self.hstr)
        if self.stdhdr[22] > 0:
            self.buf = np.fromfile(self.f, dtype=np.float32, count=(self.hdrinfo[2] - self.hdrinfo[1])//4, sep= "").byteswap()
        else:
            self.buf = np.fromfile(self.f, dtype=np.int16, count=(self.hdrinfo[2] - self.hdrinfo[1])//2, sep= "").byteswap().astype(np.float32)
        setstdparams(self)
        setspsparams(self)
        setriparams(self)
        setprparams(self)
        return 1

#  method to advance n records
    def jumpn(self, nad):
        n = nad
        while nad > 0:
            loc = self.f.tell()
            hdstr = self.f.read(12)
            self.hdrinfo = unpack_from('>4sii', hdstr)
            self.f.seek(self.hdrinfo[2] - 12, 1)
            nad -= 1

# method to go back to beginning of the file
    def rewind(self):
        self.f.seek(0, 0)

# method to move the file pointer to the next (often the first) occurrence of 'hdr_'
    def seekmark(self):
        s = self.f.read(4)
        while s != b'hdr_':
            s = self.f.read(4)
            if len(s) == 0:
                print('EOF encountered; no hdr mark')
                return 0
        self.f.seek(-4,1)
        return 1

#-------------------------------------------------------------------------------
def setspsparams(o):
    o.ipp = o.spshdr[1]
    o.gw  = o.spshdr[2]
    o.baudlen = o.spshdr[3]
    o.rf  = o.spshdr[5]
    o.codetype = o.spshdr[30]
    o.txsams = o.spshdr[35]
    o.nsamwin = o.spshdr[36]
    o.gd  = o.spshdr[37]
    o.nsam1 = o.spshdr[38]
    o.nnsam1 = o.spshdr[39]
    o.gd2 = o.spshdr[40]
    o.nsam2 = o.spshdr[41]
    o.nnsam2 = o.spshdr[42]

def setstdparams(o):
    o.date = o.stdhdr[5]
    o.time = o.stdhdr[6]
    o.hdrlen = o.stdhdr[1]
    o.reclen = o.stdhdr[2]
    o.az =  float(o.stdhdr[37])/10000.
    o.zag = float(o.stdhdr[38])/10000.
    o.zac = float(o.stdhdr[39])/10000.

def setriparams(o):
    o.bitspersam = o.rihdr[2]

def setprparams(o):
    if o.ptype == b'mracf   ' or o.ptype == b'topside ':
        o.ntimes = o.prhdr[2]
        o.nhts = o.prhdr[5]
        o.nfft = o.prhdr[6]
    elif o.ptype == b'pwr     ':
        o.ntimes = o.prhdr[2]
    elif o.ptype == b'clp     ':
        o.ntimes = o.prhdr[2]
        o.nhts = o.prhdr[4]
        o.nfft = o.prhdr[10]
        o.nadd = o.prhdr[7]
    elif o.ptype == b'dspc    ':
        o.ntimes = o.prhdr[2]
        o.nhts = o.prhdr[5]
        o.nfft = o.prhdr[6]


# Class for reading records written by the VME data taking system, waiting at EOF
# and checking every second to see if more has been written.
# This inherits from class rdvme, but overrides rdvmerec, and others.
class rdvmew(rdvme):
    def __init__(self, filestring):
        self.f = open(filestring, 'rb')
        self.filestring = filestring
        self.curloc = 0
        self.recno = 0

# method for reading a vme data taking record
# It knows about the standard data taking programs.
# (pwr, mracf, clp, topsd, dspc; these have floating point data output;
#  and rawdat; it has int16 output that is converted to float32)
# Float is read if the program hdr exists (length > 0)
# Otherwise int16 is read
# If end of file is reached, it waits for the file to grow.
    def rdvmerec(self):
        self.f.seek(self.curloc, 0)
        self.infostr = b''
        while self.infostr == b'':
            self.infostr = self.f.read(12)
            if self.infostr == b'':
                self.f.seek(self.curloc, 0)
                time.sleep(1.)
        self.hdrinfo = unpack_from('>4sii', self.infostr)
        self.f.seek(self.curloc, 0)
        self.hstr = b''
        while self.hstr == b'':
            self.hstr = self.f.read(self.hdrinfo[1])
            if self.hstr == b'':
                self.f.seek(self.curloc, 0)
                time.sleep(1.)
        self.curloc += self.hdrinfo[2]
        self.stdstr = '>4s2i8s7i16B16i'
        self.stdhdr = unpack_from(self.stdstr, self.hstr)
        self.ptype = self.stdhdr[3]
        self.ristr =  '>' + str(self.stdhdr[12]*16) + 'x' + str(self.stdhdr[14]//4) + 'i'
        self.prstr =  '>' + str(self.stdhdr[20]*16) + 'x' + str(self.stdhdr[22]//4) + 'i'
        self.spsstr = '>' + str(self.stdhdr[24]*16 + self.stdhdr[25]//16) + 'x8s7fi21f3sxi3f2if2if2i'
        self.rihdr =  unpack_from(self.ristr, self.hstr)
        cloc = self.f.tell()
        reddata = 0
        while reddata == 0:
            try:
                if self.stdhdr[22] > 0: 
                    self.prhdr = unpack_from(self.prstr, self.hstr)
                self.spshdr =  unpack_from(self.spsstr, self.hstr)
                if self.stdhdr[22] > 0:
                    self.buf = np.fromfile(self.f, dtype=np.float32, count=(self.hdrinfo[2] - self.hdrinfo[1])//4, sep= "").byteswap()
                    if (self.hdrinfo[2] - self.hdrinfo[1])//4 == self.buf.size:
                        reddata = 1
                    else:
                        self.f.seek(cloc, 0)
                        time.sleep(1.)
                else:
                    self.buf = np.fromfile(self.f, dtype=np.int16, count=(self.hdrinfo[2] - self.hdrinfo[1])//2, sep= "").byteswap().astype(np.float32)
                    if (self.hdrinfo[2] - self.hdrinfo[1])//2 == self.buf.size:
                        reddata = 1
                    else:
                        self.f.seek(cloc, 0)
                        time.sleep(1.)
            except IOError:
                self.f.seek(cloc, 0)
                time.sleep(1.)
        setstdparams(self)
        setspsparams(self)
        setriparams(self)
        setprparams(self)
        if self.recno < 2:
            print('w:', self.recno,  self.codetype, self.buf.size)
        self.recno += 1
        return 1

#  method to advance n records
    def jumpn(self, nad):
        n = nad
        while nad > 0:
            loc = self.f.tell()
            hdstr = self.f.read(12)
            self.hdrinfo = unpack_from('>4sii', hdstr)
            self.f.seek(self.hdrinfo[2] - 12, 1)
            nad -= 1
        self.curloc = self.f.tell()

# method to go back to beginning of the file
    def rewind(self):
        self.f.seek(0, 0)
        self.curloc = 0
        
# method to move the file pointer to the next (often the first) occurrence of 'hdr_'
    def seekmark(self):
        s = self.f.read(4)
        while s != b'hdr_':
            s = self.f.read(4)
            if len(s) == 0:
                print('EOF encountered; no hdr mark')
                return 0
        self.f.seek(-4,1)
        self.curloc = self.f.tell()

# Class like rdvmew, but it takes a file base name as an argument and a starting
# file number.  (The base name must include "." if needed.)
# It reads records from this file, but after it has not grown for
# a while, it checks to see if it can read from the next file in the sequence,
# continuing to tyr to read from the current file until it can access the next.
# This class inherits from rdvmew.
class rdvmewf(rdvmew):
    def __init__(self, filebase, nstart):
        self.filebase = filebase
        self.nstart = nstart
        self.ncfile = nstart
        self.filestring = filebase + '{0:03d}'.format(self.ncfile)
        self.f = open(self.filestring, 'rb')
        print(self.filestring)
        self.recno = 0
        self.curloc = 0

# method for reading a vme data taking record
# It knows about the standard data taking programs.
# (pwr, mracf, clp, topsd, dspc; these have floating point data output;
#  and rawdat; it has int16 output that is converted to float32)
# Float is read if the program hdr exists (length > 0)
# Otherwise int16 is read
# If end of file is reached, it waits for the file to grow, but
# if file does not grow after a while, the file number is incremented
# and and the routine reads the new file by calling
# rdvmew.rdvmerec  from incfileandread(self) since we want to wait
# for a record to be written in the file, but we do not want increament
# the file again if this takes a while.  That is, we wait for however
# long it takes for a record to appear in the new file.
# There are few local variables; class instance variables allow recursion.
    def rdvmerec(self):
        self.f.seek(self.curloc, 0)
        self.waitcount = 0
        self.infostr = b''
        while self.infostr == b'' or len(self.infostr) != 12:
            self.infostr = self.f.read(12)
            if self.infostr == b'' or len(self.infostr) != 12:
                self.f.seek(self.curloc, 0)
                time.sleep(1.)
                self.waitcount += 1
                if self.waitcount == 5:
                    if(incfileandread(self)):
                        return
                    else:
                        self.waitcount = 0
        self.hdrinfo = unpack_from('>4sii', self.infostr)
        self.f.seek(self.curloc, 0)
        self.waitcount = 0
        self.hstr = b''
        while self.hstr == b'' or len(self.hstr) != self.hdrinfo[1]:
            if self.hdrinfo[1] > 0 or self.hdrinfo[1] == -1:
                self.hstr = self.f.read(self.hdrinfo[1])
            if self.hstr == b'' or len(self.hstr) != self.hdrinfo[1]:
                self.f.seek(self.curloc, 0)
                time.sleep(1.)
                self.waitcount += 1
                if self.waitcount == 5:
                    if(incfileandread(self)):
                        return
                    else:
                        self.waitcount = 0
        self.curloc += self.hdrinfo[2]
        self.stdstr = '>4s2i8s7i16B16i'
        self.stdhdr = unpack_from(self.stdstr, self.hstr)
        self.ptype = self.stdhdr[3]
        self.ristr =  '>' + str(self.stdhdr[12]*16) + 'x' + str(self.stdhdr[14]//4) + 'i'
        self.prstr =  '>' + str(self.stdhdr[20]*16) + 'x' + str(self.stdhdr[22]//4) + 'i'
        self.spsstr = '>' + str(self.stdhdr[24]*16 + self.stdhdr[25]//16) + 'x8s7fi21f3sxi3f2if2if2i'
        self.rihdr =  unpack_from(self.ristr, self.hstr)
        cloc = self.f.tell()
        reddata = 0
        self.waitcount = 0
        while reddata == 0:
            try:
                if self.stdhdr[22] > 0: 
                    self.prhdr = unpack_from(self.prstr, self.hstr)
                self.spshdr =  unpack_from(self.spsstr, self.hstr)
                if self.stdhdr[22] > 0:
                    self.buf = np.fromfile(self.f, dtype=np.float32, count=(self.hdrinfo[2] - self.hdrinfo[1])//4, sep= "").byteswap()
                    if (self.hdrinfo[2] - self.hdrinfo[1])//4 == self.buf.size:
                        reddata = 1
                    else:
                        self.f.seek(cloc, 0)
                        time.sleep(1.)
                else:
                    self.buf = np.fromfile(self.f, dtype=np.int16, count=(self.hdrinfo[2] - self.hdrinfo[1])//2, sep= "").byteswap().astype(np.float32)
                    if (self.hdrinfo[2] - self.hdrinfo[1])//2 == self.buf.size:
                        reddata = 1
                    else:
                        self.f.seek(cloc, 0)
                        time.sleep(1.)
            except IOError:
                self.f.seek(cloc, 0)
                time.sleep(1.)
            self.waitcount += 1
            if self.waitcount == 5:
                if(incfileandread(self)):
                    return
                else:
                    self.waitcount = 0
        setstdparams(self)
        setspsparams(self)
        setriparams(self)
        setprparams(self)
        if self.recno < 2:
            print('wf:', self.recno,  self.codetype, self.buf.size)
        self.recno += 1
        return 1

def incfileandread(obj):
    try:
        filestring = obj.filebase + '{0:03d}'.format(obj.ncfile + 1)
        tempf = open(filestring, 'rb')
        tempf.close()
    except IOError:
        return 0
    obj.ncfile += 1
    obj.filestring = obj.filebase + '{0:03d}'.format(obj.ncfile)
    obj.f.close()
    obj.f = open(obj.filestring, 'rb')
    obj.curloc = 0
    print('New file:', obj.filestring)
    obj.recno = 0
    obj.seekmark()
    obj.rdvmerec()

# Each record of the input file is read using an instance of rdvme. 
# Then, if inf.ptype is equa to the third argument, it is
# read again as a string for convenient writing to the output file
def copyptype(infilestr, outfilestr, ptype): 
        inf = rdvme(infilestr)
        inf.rewind()
        outfile = open(outfilestr, 'wb', 0)
        readrec = 1
        while readrec:
            readrec = inf.rdvmerec()
            if readrec:
                if ptype == inf.ptype:
                    inf.f.seek(-inf.reclen, 1)
                    bbuf = inf.f.read(inf.reclen)
                    outfile.write(bbuf)
            else:
                inf.f.close()
                outfile.close()


# For testing:

#  Function copyslow opens infilestr, reads records every
# three seconds, writing out each to outfilestr just after it is read.
# For example: rdvme.copyslow(inf, of)
# makes a growing file. If  the of is opened in another window as df
# with class rdvmew then the statement:
# res = 1
# while res: res = df.rdvmerec(); print df.time, df.reclen, amax(df.buf)
# results in reading the records as they are written.
# Each record of the in file is read using an instance of rdvme.  Then it is
# read again as a string for convenient writing to the output file
def copyslow(infilestr, outfilestr): 
        inf = rdvme(infilestr)
        outfile = open(outfilestr, 'wb', 0)
        readrec = 1
        print(readrec)
        while readrec:
            readrec = inf.rdvmerec()
            print(inf.time, inf.reclen, np.amax(inf.buf))
            if readrec:
                inf.f.seek(-inf.reclen, 1)
                bbuf = inf.f.read(inf.reclen)
                outfile.write(bbuf)
                time.sleep(3.)
            else:
                inf.f.close()
                outfile.close()

# like cpoyslow, but makes multiple output files
def copymf(infilestr, outbase): 
        inf = rdvme(infilestr)
        nfile = 0
        outfilestr = outbase  + '{0:03d}'.format(nfile)
        outfile = open(outfilestr, 'wb', 0)
        nrecsperfile = 10
        nreccfile = 0
        readrec = 1
        print(readrec)
        while readrec:
            readrec = inf.rdvmerec()
            print(inf.time, inf.reclen, np.amax(inf.buf))
            if readrec:
                inf.f.seek(-inf.reclen, 1)
                bbuf = inf.f.read(inf.reclen)
                outfile.write(bbuf)
                nreccfile += 1
                if nreccfile == nrecsperfile:
                    outfile.close()
                    nfile += 1
                    outfilestr = outbase  + '{0:03d}'.format(nfile)
                    outfile = open(outfilestr, 'wb')
                    nreccfile = 0
                time.sleep(3.0)
            else:
                inf.f.close()
                outfile.close()

# copy Some records from input to output files
def copyn(infilestr, outfilestr, n): 
        inf = rdvme(infilestr)
        outfile = open(outfilestr, 'wb')
        readrec = 1
        print(readrec)
        nrec = n
        while nrec > 0:
            readrec = inf.rdvmerec()
            print(inf.time, inf.reclen, np.amax(inf.buf))
            if readrec:
                inf.f.seek(-inf.reclen, 1)
                bbuf = inf.f.read(inf.reclen)
                outfile.write(bbuf)
            else:
                inf.f.close()
                outfile.close()
            nrec -= 1
        inf.f.close()
        outfile.close()
