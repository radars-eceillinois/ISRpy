#
#  kmltools.py
#
#	tools for generating kml files for plotting over google Earth
#
#  Created by Pablo Reyes on 08/30/13 as kmltools.py
#  Copyright (c) 2013 ECE, UIUC. All rights reserved.
#  history
#  - Aug30,2013 by P. Reyes
#    -File kmltools created with class for creating a file
#  - June 5, 2016 by E. Kudeki
#    -importing radarbeam from beampack rather than radarpack
#    -beampack is now a sub package of ISRpy.


# ------------ radar specifications -------------------------
import numpy as np
import StringIO as __StringIO__

class kmzfile:
    """Will contain functions to create a kmz file"""
    def __init__(self, fname=None, title = None):
        import zipfile as __zipfile__
        self.fname = fname
        self.title = title
        self.main_kml = kmlfile('doc.kml',title = title)
        self.tempfile = __StringIO__.StringIO()
        self.zf = __zipfile__.ZipFile(self.tempfile,'w',__zipfile__.ZIP_DEFLATED)

    def __closefile__(self):
        self.zf.close()

    def savefile(self):
        self.__closefile__()
        self.tempfile.seek(0)
        fp = open(self.fname,'w')
        fp.write(self.tempfile.read())
        fp.close()
        self.tempfile.close()

    def __get_model_path__(self):
        if not hasattr(self, "__model_count__"):
            self.__model_count__ = 0
        return 'model%d.dae'%self.__model_count__

    def __get_image_path__(self):
        if not hasattr(self, "__image_count__"):
            self.__image_count__ = 0
        return 'image%d.png'%self.__image_count__

    def add_rectangular_image(self, image, width, height, lon, lat, ht,
        heading =0.0, tilt=0.0, roll=0.0, rad=False,transparency=0.1):
        self.main_kml.lookat(lon, lat, ht, tilt=np.pi/4,
                             range_dist=5*max(width,height)/1000.,rad=rad)
        dae_path = self.__get_model_path__()
        self.main_kml.placemodel(lon, lat, ht, dae_path=dae_path, rad=rad,
                                 heading=heading, tilt=tilt, roll=roll)
        self.main_kml.__close_kml_file__()
        self.zf.writestr(zinfo_or_arcname='doc.kml', bytes=self.main_kml.outstr)
        fp_daemodel = __StringIO__.StringIO() # starting file instance
        img_path = self.__get_image_path__() # getting current image path name
        rect_vertical_model(fp_daemodel,img_path,width=width,height=height,
                            transparency=transparency)
        fp_daemodel.seek(0) # ready to be read
        self.zf.writestr(zinfo_or_arcname=dae_path, bytes= fp_daemodel.read())
        fp_daemodel.close()
        if image.__class__.__name__ == 'Figure':
            fp_image = __StringIO__.StringIO()
            image.savefig(fp_image)
            fp_image.seek(0)
            self.zf.writestr(zinfo_or_arcname=img_path, bytes= fp_image.read())
            fp_image.close()
        elif image.__class__.__name__ == 'str':
            self.zf.write(image, arcname=img_path)
            pass

    def add_curved_image(self, image, lon, lat, ht, r0=100000., r1=200000.,
        theta = np.arange(0.,31)*np.pi/90, phi = np.ones(31)* 0.,
        heading =0.0, tilt=0.0, roll=0.0, rad=False,transparency=0.1):
        """
        theta = np.arange(0.,31)*np.pi/90
        phi = np.ones(31)* 0.
        """
        self.main_kml.lookat(lon, lat, ht, tilt=np.pi/4,
                             range_dist=5*max(r0,r1)/1000.,rad=True)
        dae_path = self.__get_model_path__()
        self.main_kml.placemodel(lon, lat, ht, dae_path=dae_path, rad=rad,
                                 heading=heading, tilt=tilt, roll=roll)
        self.main_kml.__close_kml_file__()
        self.zf.writestr(zinfo_or_arcname='doc.kml', bytes=self.main_kml.outstr)
        fp_daemodel = __StringIO__.StringIO() # starting file instance
        img_path = self.__get_image_path__() # getting current image path name
        radial_surface_model(fp_daemodel,img_path,r0=r0,r1=r1,
                             theta=theta, phi=phi,
                             transparency=transparency)
        fp_daemodel.seek(0) # ready to be read
        self.zf.writestr(zinfo_or_arcname=dae_path, bytes= fp_daemodel.read())
        fp_daemodel.close()
        if image.__class__.__name__ == 'Figure':
            fp_image = __StringIO__.StringIO()
            image.savefig(fp_image)
            fp_image.seek(0)
            self.zf.writestr(zinfo_or_arcname=img_path, bytes= fp_image.read())
            fp_image.close()
        elif image.__class__.__name__ == 'str':
            self.zf.write(image, arcname=img_path)
            pass


class kmlfile:
    """Will contain functions to create a kml file"""
    def __init__(self,fname=None,title=None):
        self.fname = fname
        self.outstr =  """<?xml version="1.0" encoding="UTF-8"?>\n"""
        self.outstr += """<kml xmlns="http://www.opengis.net/kml/2.2"\n"""
        self.outstr += """xmlns:gx="http://www.google.com/kml/ext/2.2"\n"""
        self.outstr += """xmlns:kml="http://www.opengis.net/kml/2.2"\n"""
        self.outstr += """xmlns:atom="http://www.w3.org/2005/Atom">\n"""
        self.outstr +=   "<Document>\n"
        if title!=None:
            self.outstr += "<name>"+title+"</name>\n"

    def __close_kml_file__(self):
        self.outstr += "</Document>\n"
        self.outstr += "</kml>\n"

    def savefile(self, fname=None):
        self.__close_kml_file__()
        if fname != None:
            self.fname = fname
        if self.fname==None:
            import os
            basen = 'Untitled'
            i=0
            while os.path.exists(basen+str(i)+'.kml'):
                i += 1
            self.fname = basen+str(i)+'.kml'
        if type(self.fname) == str:
            try:
                fp = open(self.fname,'w')
                fp.write(self.outstr)
                fp.close()
            except:
                print "Error writing  ", self.fname
        else:
            try:
                self.fname.write(self.outstr)
            except:
                print "Error writing in ", self.fname

    def startfolder(self,foldertitle=None):
        self.outstr += "<Folder>\n"
        if foldertitle!=None:
            self.outstr += "<name>"+foldertitle+"</name>\n"

    def addtime(self,yyyy,mm,dd,HH=0,MM=0,SS=0):
        self.outstr += "<TimeStamp>\n"
        self.outstr += "  <when>"
        self.outstr += "%.4d-%.2d-%.2dT%.2d:%.2d:%.2d"%(yyyy,mm,dd,HH,MM,SS)
        self.outstr += "</when>\n"
        self.outstr += "</TimeStamp>\n"

    def addtimespan(self,yyyy1,mm1,dd1,yyyy2,mm2,dd2,
            HH1=0,MM1=0,SS1=0,HH2=23,MM2=59,SS2=59):
        self.outstr += "<TimeSpan>\n"
        self.outstr += "  <begin>"
        self.outstr += "%.4d-%.2d-%.2dT%.2d:%.2d:%.2d"%(yyyy1,mm1,dd1,HH1,MM1,SS1)
        self.outstr += "</begin>\n"
        self.outstr += "  <end>"
        self.outstr += "%.4d-%.2d-%.2dT%.2d:%.2d:%.2d"%(yyyy2,mm2,dd2,HH2,MM2,SS2)
        self.outstr += "</end>\n"
        self.outstr += "</TimeSpan>\n"

    def closefolder(self):
        self.outstr += "</Folder>\n"

    def plotline(self,lons,lats,hts,linestyle=None,color=None,width=None,
                 name=None,rad=False,extrude=0,yyyy=None,mm=None,dd=None,
                 yyyy2=None,mm2=None,dd2=None):
        """ plotline(lons,lats,hts,linestyle=None,color=None,width=None,
                name=None,rad=False,extrude=0)
            degs, degs, km
            color: alpha,B,G,R, example "ff0000ff":red
            width: default=1
            name: Name for the collection of points that makes the line
            rad: if True the lons and lats must be in radians,
                 degrees are used otherwise
            extrude: if 1 a line is draw to the ground
        """
        points = len(lons)
        outp = "<Placemark>\n"
        if name != None:
            outp += "  <name>%s</name>\n"%name
        if yyyy!=None and mm!=None and dd!=None:
            if yyyy2!=None and mm2!=None and dd2!=None:
                outp += "<TimeSpan>\n"
                outp += "  <begin>%.4d-%.2d-%.2d</begin>\n"%(yyyy,mm,dd)
                outp += "  <end>%.4d-%.2d-%.2d</end>\n"%(yyyy2,mm2,dd2)
                outp += "</TimeSpan>\n"
            else:
                outp += "<TimeStamp>\n"
                outp += "  <when>%.4d-%.2d-%.2d</when>\n"%(yyyy,mm,dd)
                outp += "</TimeStamp>\n"

        if linestyle!=None:
            outp += """  <styleUrl>%s</styleUrl>\n"""%linestyle
        else:
            outp += "  <Style>\n"
            outp += "    <LineStyle>\n"
            if color!=None:
                outp += "      <color>%s</color>\n"%color
            if width!=None:
                outp += "      <width>%s</width>\n"%width
            outp += "    </LineStyle>\n"
            outp += "  </Style>\n"
        outp += "  <LineString>\n"
        outp += "    <extrude>%d</extrude>\n"%extrude
        outp += "    <altitudeMode>absolute</altitudeMode>\n"
        outp += "    <coordinates>\n"
        outp += "      "
        for  i in range(points):
            if rad:
                outp += "%f,%f,%f "%(lons[i]*180./np.pi,lats[i]*180./np.pi,1000 * hts[i])
            else:
                outp += "%f,%f,%f "%(lons[i],lats[i],1000 * hts[i])
        outp += "\n"
        outp += "    </coordinates>\n"
        outp += "  </LineString>\n"
        outp += "</Placemark>\n"
        self.outstr += outp

    def lookat(self, lon,lat,ht,heading=0.0, tilt=0.0, range_dist = 1.0, rad=False):
        if rad:
            lon0,lat0,alt0 = (lon*180./np.pi, lat*180./np.pi, 1000 * ht)
            heading0,tilt0 = (heading*180./np.pi, tilt*180./np.pi)
        else:
            lon0,lat0,alt0 = (lon,lat,1000 * ht)
            heading0,tilt0 = (heading, tilt)
        self.outstr += "<LookAt id=\"ID0\">\n"
        self.outstr += "  <longitude>%f</longitude>\n"%lon0
        self.outstr += "  <latitude>%f</latitude>\n"%lat0
        self.outstr += "  <altitude>%f</altitude>\n"%alt0
        self.outstr += "  <heading>%f</heading>\n"%heading0
        self.outstr += "  <tilt>%f</tilt>\n"%tilt0
        self.outstr += "  <range>%f</range>\n"%(range_dist * 1000.0)
        self.outstr += "  <altitudeMode>relativeToGround</altitudeMode>\n"
        self.outstr += "</LookAt>\n"

    def camera(self, lon,lat,ht,heading=0.0, tilt=0.0, roll = 0.0, rad=False):
        if rad:
            lon0,lat0,alt0 = (lon*180./np.pi, lat*180./np.pi, 1000 * ht)
            heading0,tilt0,roll0 = (heading*180./np.pi, tilt*180./np.pi, roll*180./np.pi)
        else:
            lon0,lat0,alt0 = (lon,lat,1000 * ht)
            heading0,tilt0,roll0 = (heading, tilt, roll)
        self.outstr += "<Camera id=\"ID1\">\n"
        self.outstr += "  <longitude>%f</longitude>\n"%lon0
        self.outstr += "  <latitude>%f</latitude>\n"%lat0
        self.outstr += "  <altitude>%f</altitude>\n"%alt0
        self.outstr += "  <heading>%f</heading>\n"%heading0
        self.outstr += "  <tilt>%f</tilt>\n"%tilt0
        self.outstr += "  <roll>%f</roll>\n"%roll0
        self.outstr += "  <altitudeMode>relativeToGround</altitudeMode>\n"
        self.outstr += "</Camera>\n"

    def placemodel(self,lon,lat,ht, dae_path, name=None, rad=False,
            heading=0.0, tilt=0.0, roll=0.0, scalexyz=(1.0,1.0,1.0),
            resources=['2Dimage.png']):
        outp = "<Placemark>\n"
        if name != None:
            outp += "  <name>%s</name>\n"%name
        outp += "  <Model id=\"ID\">\n"
        outp += "    <altitudeMode>relativeToGround</altitudeMode>\n"
        outp += "    <Location>\n"
        if rad:
            lon0,lat0,alt0 = (lon*180./np.pi, lat*180./np.pi, 1000 * ht)
        else:
            lon0,lat0,alt0 = (lon,lat,1000 * ht)
        outp += "      <longitude>%f</longitude>\n"%lon0
        outp += "      <latitude>%f</latitude>\n"%lat0
        outp += "      <altitude>%f</altitude>\n"%alt0
        outp += "    </Location>\n"
        outp += "    <Orientation>\n"
        outp += "      <heading>%f</heading>\n"%heading
        outp += "      <tilt>%f</tilt>\n"%tilt
        outp += "      <roll>%f</roll>\n"%roll
        outp += "    </Orientation>\n"
        outp += "    <Scale>\n"
        outp += "      <x>%f</x>\n"%scalexyz[0]
        outp += "      <y>%f</y>\n"%scalexyz[1]
        outp += "      <z>%f</z>\n"%scalexyz[2]
        outp += "    </Scale>\n"
        outp += "    <Link>\n"
        outp += "      <href>%s</href>\n"%dae_path
        outp += "    </Link>\n"
        outp += "    <ResourceMap>\n"
        for rname in resources:
            outp += "      <Alias>\n"
            outp += "      <targetHref>%s</targetHref>\n"%rname
            outp += "      <sourceHref>%s</sourceHref>\n"%rname
            outp += "      </Alias>\n"
        outp += "    </ResourceMap>\n"
        outp += "  </Model>\n"
        outp += "</Placemark>\n"
        self.outstr += outp


    def plotpoint(self,lon,lat,ht,iconstyle=None,color=None,scale=None,
                 name=None,rad=False,iconref=None,extrude=0):
        """ plotpoint(lon,lat,ht,iconstyle=None,color=None,scale=None,name=None,
                       rad=False,iconref=None,extrude=0)
            degs, degs, km
            iconstyle: used if a predetermine set of styles has been defined
            color: alpha,B,G,R, example "ff0000ff":red
            scale: default=1
            name: Name for the collection of points that makes the line
            rad: if True the lons and lats are given in radians, degrees are used otherwise
            iconref: if None, the default is placed. To check the available href icons
            right click on a placemark in google earth, select "Get Info" and click on the
            icon simbol to check the url location.
            example= "http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png"
            extrude: if 1 a line is draw to the ground
        """
        outp = "<Placemark>\n"
        if name != None:
            outp += "  <name>%s</name>\n"%name
        if iconstyle!=None:
            outp += """  <styleUrl>%s</styleUrl>\n"""%linestyle
        else:
            outp += "  <Style>\n"
            outp += "    <IconStyle>\n"
            if color!=None:
                outp += "      <color>%s</color>\n"%color
            if scale!=None:
                outp += "      <scale>%s</scale>\n"%scale
            if iconref!=None:
                outp += "      <Icon>\n"
                outp += "        <href>%s</href>\n"%iconref
                outp += "      </Icon>\n"
            outp += """      <hotSpot x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>\n"""
            outp += "    </IconStyle>\n"
            outp += "  </Style>\n"
        outp += "  <Point>\n"
        outp += "    <extrude>%d</extrude>\n"%extrude
        outp += "    <altitudeMode>absolute</altitudeMode>\n"
        outp += "    <coordinates>\n"
        outp += "      "
        if rad:
            outp += "%f,%f,%f "%(lon*180./np.pi,lat*180./np.pi,1000 * ht)
        else:
            outp += "%f,%f,%f "%(lon,lat,1000 * ht)
        outp += "\n"
        outp += "    </coordinates>\n"
        outp += "  </Point>\n"
        outp += "</Placemark>\n"
        self.outstr += outp

    def plotPolygon(self,outlons,outlats,outhts,linecolor=None,linewidth=None,
                    fillcolor=None, name=None,Style=None,
                    inlons=None,inlats=None,inhts=None,extrude=0,rad=False):
        """
        plotcirperp2V(outlons,outlats,outhts,linecolor=None,linewidth=None,
                    fillcolor=None, name=None,Style=None,
                    inlons=None,inlats=None,inhts=None,extrude=0,rad=False)
        """
        outp = "<Placemark>\n"
        if name != None:
            outp += "  <name>%s</name>\n"%name
        if Style!=None:
            outp += """  <styleUrl>%s</styleUrl>\n"""%Style
        else:
            outp += "  <Style>\n"
            if linecolor!=None:
                outp += "    <LineStyle>\n"
                outp += "      <color>%s</color>\n"%linecolor
                if linewidth!=None:
                    outp += "      <width>%s</width>\n"%linewidth
                outp += "    </LineStyle>\n"
            if fillcolor!=None:
                outp += "    <PolyStyle>\n"
                outp += "      <color>%s</color>\n"%fillcolor
                outp += "    </PolyStyle>\n"
            outp += "  </Style>\n"
        outp += "  <Polygon>\n"
        outp += "    <extrude>%d</extrude>\n"%extrude
        outp += "    <altitudeMode>absolute</altitudeMode>\n"
        outp += "    <outerBoundaryIs>\n"
        outp += "     <LinearRing>\n"
        outp += "      <coordinates>\n"
        outp += "        "
        for  i in range(len(outlons)):
            if rad:
                outp += "%f,%f,%f "%(outlons[i]*180./np.pi,outlats[i]*180./np.pi,1000. * outhts[i])
            else:
                outp += "%f,%f,%f "%(outlons[i],outlats[i],1000. * outhts[i])
        outp += "\n"
        outp += "       </coordinates>\n"
        outp += "      </LinearRing>\n"
        outp += "    </outerBoundaryIs>\n"
        if inlons!=None and inlats!=None and inhts!=None:
            outp += "    <innerBoundaryIs>\n"
            outp += "     <LinearRing>\n"
            outp += "      <extrude>%d</extrude>\n"%extrude
            outp += "      <altitudeMode>absolute</altitudeMode>\n"
            outp += "      <coordinates>\n"
            outp += "        "
            for  i in range(len(inlats)):
                if rad:
                    outp += "%f,%f,%f "%(inlons[i]*180./np.pi,inlats[i]*180./np.pi,1000. * inhts[i])
                else:
                    outp += "%f,%f,%f "%(inlons[i],inlats[i],1000. * inhts[i])
            outp += "\n"
            outp += "       </coordinates>\n"
            outp += "      </LinearRing>\n"
            outp += "    </innerBoundaryIs>\n"

        outp += "  </Polygon>\n"
        outp += "</Placemark>\n"
        self.outstr += outp

    def plot_square_perp2V(self,lon,lat,htkm,Veast, Vnorth, Vup,
                    Vx=None,Vy=None,Vz=None,sidelen=1.,
                    linecolor=None,linewidth=None,
                    fillcolor=None, name=None,Style=None,extrude=0,rad=False):
        """
        plot_square_perp2V(lon,lat,htkm,Veast, Vnorth, Vup,sidelen=1.,
                    linecolor=None,linewidth=None,
                    fillcolor=None, name=None,Style=None,extrude=0,rad=False)
        """
        from beampack import radarbeam as rb
        Xr,Yr,Zr = rb.llh2xyz(lat, lon, htkm) # current point
        if Vx==None and Vy==None and Vz==None:
            X,Y,Z = rb.enu2xyz(Xr,Yr,Zr,Veast,Vnorth,Vup)
            w = np.array([X-Xr,Y-Yr,Z-Zr])
        else:
            w = np.array([Vx,Vy,Vz])

        # finding a couple of orthogonal unit vectors perpendicular to w
        # such that w/|w| = u x v
        t = np.zeros_like(w)
        maxind = np.argmax(w)
        t[(maxind+1) % w.size ] = w[maxind]
        u = np.cross(w,t)
        u = u/ np.linalg.norm(u)
        v = np.cross(w,u)
        v = v/ np.linalg.norm(v)

        X,Y,Z = np.array([Xr,Yr,Zr]) + sidelen/2.* u + sidelen/2.* v
        la1,lo1,ht1=rb.xyz2llh(X,Y,Z)
        X,Y,Z = np.array([Xr,Yr,Zr]) - sidelen/2.* u + sidelen/2.* v
        la2,lo2,ht2=rb.xyz2llh(X,Y,Z)
        X,Y,Z = np.array([Xr,Yr,Zr]) - sidelen/2.* u - sidelen/2.* v
        la3,lo3,ht3=rb.xyz2llh(X,Y,Z)
        X,Y,Z = np.array([Xr,Yr,Zr]) + sidelen/2.* u - sidelen/2.* v
        la4,lo4,ht4=rb.xyz2llh(X,Y,Z)
        outlons = [lo1,lo2,lo3,lo4,lo1]
        outlats = [la1,la2,la3,la4,la1]
        outhts = [ht1,ht2,ht3,ht4,ht1]
        self.plotPolygon(outlons, outlats, outhts, linecolor=linecolor,
            linewidth=linewidth, fillcolor=fillcolor, name=name, Style=Style,
            inlons=None, inlats=None, inhts=None, extrude=extrude, rad=rad)

    def plot_circle_perp2V(self, lon, lat, htkm, Veast=None, Vnorth=None, Vup=None,
            Vx=None,Vy=None,Vz=None,
            diameter=1., linecolor=None, linewidth=None, fillcolor=None,
            name=None, Style=None, extrude=0, rad=False):
        """
        plot_square_perp2V(lon,lat,htkm,Veast, Vnorth, Vup,diameter=1.,
                    linecolor=None,linewidth=None,
                    fillcolor=None, name=None,Style=None,extrude=0,rad=False)
        """
        from beampack import radarbeam as rb
        Xr,Yr,Zr = rb.llh2xyz(lat, lon, htkm) # current point
        if Vx==None and Vy==None and Vz==None:
            X,Y,Z = rb.enu2xyz(Xr,Yr,Zr,Veast,Vnorth,Vup)
            w = np.array([X-Xr,Y-Yr,Z-Zr])
        else:
            w = np.array([Vx,Vy,Vz])
        # finding a couple of orthogonal unit vectors perpendicular to w
        # such that w/|w| = u x v
        t = np.zeros_like(w)
        maxind = np.argmax(w)
        t[(maxind+1) % w.size ] = w[maxind]
        u = np.cross(w,t)
        u = u/ np.linalg.norm(u)
        v = np.cross(w,u)
        v = v/ np.linalg.norm(v)

        outlons = []
        outlats = []
        outhts = []
        r = diameter/2.
        for phi in np.arange(0,2*np.pi+np.pi/10.,np.pi/10.):
            X,Y,Z = np.array([Xr,Yr,Zr]) + u*r*np.cos(phi) + v*r*np.sin(phi)
            lat1,lon1,ht1=rb.xyz2llh(X,Y,Z)
            outlons += [lon1]
            outlats += [lat1]
            outhts += [ht1]

        self.plotPolygon(outlons,outlats,outhts,linecolor=linecolor,
            linewidth=linewidth, fillcolor=fillcolor, name=name, Style=Style,
            inlons=None, inlats=None, inhts=None, extrude=extrude, rad=rad)

        """
        X,Y,Z = np.array([Xr,Yr,Zr]) + sidelen/2.* u + sidelen/2.* v
        la1,lo1,ht1=rb.xyz2llh(X,Y,Z)
        X,Y,Z = np.array([Xr,Yr,Zr]) - sidelen/2.* u + sidelen/2.* v
        la2,lo2,ht2=rb.xyz2llh(X,Y,Z)
        X,Y,Z = np.array([Xr,Yr,Zr]) - sidelen/2.* u - sidelen/2.* v
        la3,lo3,ht3=rb.xyz2llh(X,Y,Z)
        X,Y,Z = np.array([Xr,Yr,Zr]) + sidelen/2.* u - sidelen/2.* v
        la4,lo4,ht4=rb.xyz2llh(X,Y,Z)
        outlons = [lo1,lo2,lo3,lo4,lo1]
        outlats = [la1,la2,la3,la4,la1]
        outhts = [ht1,ht2,ht3,ht4,ht1]
        self.plotPolygon(outlons,outlats,outhts,linecolor=linecolor,linewidth=linewidth,
                    fillcolor=fillcolor, name=name,Style=Style,
                    inlons=None,inlats=None,inhts=None,extrude=extrude,rad=rad)
        """

def rect_vertical_model(fname, imagename, width = 1000, height = 1000, transparency=0.1):
    import collada as cda
    mesh = cda.Collada()
    axis = cda.asset.UP_AXIS.Z_UP
    mesh.assetInfo.upaxis = axis

    image1 = cda.material.CImage("material_img", imagename)
    surface = cda.material.Surface("material-surface", image1)
    sampler2d = cda.material.Sampler2D("material-sampler",surface)
    map1 = cda.material.Map(sampler2d, "sampler_1")

    effect1 = cda.material.Effect("effect0", [surface, sampler2d], "lambert",
        emission=None, ambient=None, diffuse=map1,
        specular=None, shininess=None, reflective=None,
        reflectivity=None, transparent=None, transparency=transparency,
        double_sided=True)
    mat1 = cda.material.Material( id="materialID", name="2Dimage", effect=effect1)
    mesh.effects.append(effect1)
    mesh.materials.append(mat1)
    mesh.images.append(image1)

    import numpy as np

    vert_f = np.array([0, 0, 0, 0, 0, height,width, 0, height, width, 0, 0])
                       # x0, y0, z0,    x1, y1, z1,    ...
    normal_f = np.array([   0,   1,  0,      0,   -1,  0])
                        # ux0, uy0, uz0,    ux1, uy1, uz1,    ...
    textcoord_f = np.array([0,  0,     0,  1,    1,  1,    1,  0])
                        # s0, t0,    s1, t1,    ...
    ver_src = cda.source.FloatSource("verts", vert_f, ('X', 'Y', 'Z'))
    nor_src = cda.source.FloatSource("normals", normal_f, ('X', 'Y', 'Z'))
    tex_src = cda.source.FloatSource("texcoords", textcoord_f, ('S', 'T'))
    geom = cda.geometry.Geometry(mesh, "geo0", "geo0", [ver_src, nor_src, tex_src])
    input_list = cda.source.InputList()
    input_list.addInput(0, 'VERTEX', "#verts")
    input_list.addInput(1, 'NORMAL', "#normals")
    input_list.addInput(2, 'TEXCOORD', "#texcoords")
    indices = np.array([0, 0, 0,    1, 0, 1,     2, 0, 2,
                        0, 0, 0,    2, 0, 2,     3, 0, 3,
                        0, 1, 0,    2, 1, 2,     1, 1, 1,
                        0, 1, 0,    3, 1, 3,     2, 1, 2 ])
    triset = geom.createTriangleSet(indices, input_list, "material")
    geom.primitives.append(triset)
    mesh.geometries.append(geom)
    matnode = cda.scene.MaterialNode("material", mat1, inputs=[])
    geomnode = cda.scene.GeometryNode(geom, [matnode])
    node = cda.scene.Node("node0", children=[geomnode])
    myscene = cda.scene.Scene("myscene", [node])
    mesh.scenes.append(myscene)
    mesh.scene = myscene
    mesh.write(fname)

def radial_surface_model(fname, imagename, transparency=0.1, r0=100,r1=200,
            theta = np.arange(0.,31)*np.pi/90, phi = np.ones(31)* 0.):
    import collada as cda
    mesh = cda.Collada()
    axis = cda.asset.UP_AXIS.Z_UP
    mesh.assetInfo.upaxis = axis

    image1 = cda.material.CImage("material_img", imagename)
    surface = cda.material.Surface("material-surface", image1)
    sampler2d = cda.material.Sampler2D("material-sampler",surface)
    map1 = cda.material.Map(sampler2d, "sampler_1")

    effect1 = cda.material.Effect("effect0", [surface, sampler2d], "lambert",
        emission=None, ambient=None, diffuse=map1,
        specular=None, shininess=None, reflective=None,
        reflectivity=None, transparent=None, transparency=transparency,
        double_sided=True)
    mat1 = cda.material.Material( id="materialID", name="2Dimage", effect=effect1)
    mesh.effects.append(effect1)
    mesh.materials.append(mat1)
    mesh.images.append(image1)

    import numpy as np

 #   r0 = 100
 #   r1 = 200
 #   ndiv = 30 # 3 divisions
#    st_delta = 1./ndiv
    nradii = len(theta)
    ndiv = nradii - 1
#    nradii = ndiv+1 # number of vertices
#    phi = np.ones(nradii)* 0.
 #   del_theta = np.pi/90
  #  theta = np.arange(0.,nradii)*del_theta
    ntri = 4 * ndiv
    vert_f = []
    for i in range(nradii):
        u  = np.array([np.sin(theta[i])*np.cos(phi[i]),
                       np.sin(theta[i])*np.sin(phi[i]),
                       np.cos(theta[i])])
        vert_f += list(u*r0)
        vert_f += list(u*r1)
    vert_f = np.array(vert_f)
    textcoord_f = np.zeros(2*2*nradii)
    textcoord_f[3::4]=1
    textcoord_f[2::4]=textcoord_f[0::4]=np.arange(0,ndiv+1.)/ndiv
    v1s = [1,3,3,2]
    v2s = [3,1,2,3]
    normal_f = []
    indices = []
    for i in range(ntri):
        vi0 = i/4*2
        vi1 = vi0+v1s[i%4]
        vi2 = vi0+v2s[i%4]
        V0 = np.array(vert_f[3*vi0:3*(vi0+1)])
        V1 = np.array(vert_f[3*vi1:3*(vi1+1)])
        V2 = np.array(vert_f[3*vi2:3*(vi2+1)])
        n0 = np.cross(V1-V0,V2-V0)
        n0 = n0 / np.sqrt(np.sum(n0**2))
        normal_f += list(n0)
        indices += [vi0,i,vi0, vi1,i,vi1, vi2,i,vi2]
    normal_f = np.array(normal_f)
    indices = np.array(indices)


    """
    vert_f = np.array([0, 0, 0, 0, 0, height,width, 0, height, width, 0, 0])
                       # x0, y0, z0,    x1, y1, z1,    ...
    normal_f = np.array([   0,   1,  0,      0,   -1,  0])
                        # ux0, uy0, uz0,    ux1, uy1, uz1,    ...
    textcoord_f = np.array([0,  0,     0,  1,    1,  1,    1,  0])
                        # s0, t0,    s1, t1,    ...
    """
    ver_src = cda.source.FloatSource("verts", vert_f, ('X', 'Y', 'Z'))
    nor_src = cda.source.FloatSource("normals", normal_f, ('X', 'Y', 'Z'))
    tex_src = cda.source.FloatSource("texcoords", textcoord_f, ('S', 'T'))
    geom = cda.geometry.Geometry(mesh, "geo0", "geo0", [ver_src, nor_src, tex_src])
    input_list = cda.source.InputList()
    input_list.addInput(0, 'VERTEX', "#verts")
    input_list.addInput(1, 'NORMAL', "#normals")
    input_list.addInput(2, 'TEXCOORD', "#texcoords")
    """
    indices = np.array([0, 0, 0,    1, 0, 1,     2, 0, 2,
                        0, 0, 0,    2, 0, 2,     3, 0, 3,
                        0, 1, 0,    2, 1, 2,     1, 1, 1,
                        0, 1, 0,    3, 1, 3,     2, 1, 2 ])
    """
    triset = geom.createTriangleSet(indices, input_list, "material")
    geom.primitives.append(triset)
    mesh.geometries.append(geom)
    matnode = cda.scene.MaterialNode("material", mat1, inputs=[])
    geomnode = cda.scene.GeometryNode(geom, [matnode])
    node = cda.scene.Node("node0", children=[geomnode])
    myscene = cda.scene.Scene("myscene", [node])
    mesh.scenes.append(myscene)
    mesh.scene = myscene
    mesh.write(fname)
