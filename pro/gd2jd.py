#!/usr/bin/python
import sys
def gd2jd(st):
     if len(st) != 8:
         return '00000'
     yyyy, mm, dd = int(st[:4]), int(st[4:6]), int(st[6:])
     hh, min, sec = 12.0, 0.0, 0.0
 
     UT=hh+min/60+sec/3600
 
     #print "UT="+`UT`
 
     total_seconds=hh*3600+min*60+sec
     fracday=total_seconds/86400
 
     #print "Fractional day: %f" % fracday
     # print dd,mm,yyyy, hh,min,sec, UT
 
     if (100*yyyy+mm-190002.5)>0:
         sig=1
     else:
         sig=-1
 
     JD = 367*yyyy - int(7*(yyyy+int((mm+9)/12))/4) + int(275*mm/9) + dd + 1721013.5 + UT/24 - 0.5*sig +0.5
     return str(int(JD))[2:]
#print gd2jd(sys.argv[1])
