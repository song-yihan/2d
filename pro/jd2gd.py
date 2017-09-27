#!/usr/bin/python
import sys
def jd2gd(dd):
    dd = str(dd)
    #2450517

    if dd[:2] !='24':
        dd = '24'+dd

    x = dd
    try:
        jd=float(x)
    except ValueError:
        return ''

    jd=jd+0.5
    Z=int(jd)
    F=jd-Z
    alpha=int((Z-1867216.25)/36524.25)
    A=Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int( (B-122.1)/365.25)
    D = int( 365.25*C )
    E = int( (B-D)/30.6001 )

    dd = B - D - int(30.6001*E) + F

    if E<13.5:
        mm=E-1

    if E>13.5:
        mm=E-13

    if mm>2.5:
        yyyy=C-4716

    if mm<2.5:
        yyyy=C-4715

    daylist=[31,28,31,30,31,30,31,31,30,31,30,31]
    daylist2=[31,29,31,30,31,30,31,31,30,31,30,31]

    h=int((dd-int(dd))*24)
    min=int((((dd-int(dd))*24)-h)*60)
    sec=86400*(dd-int(dd))-h*3600-min*60

    # Now calculate the fractional year. Do we have a leap year?
    if (yyyy%4 != 0):
        days=daylist2
    elif (yyyy%400 == 0):
        days=daylist2
    elif (yyyy%100 == 0):
        days=daylist
    else:
        days=daylist2


    ret = '%4d%2d%2dA---A%2d:%2d:%4.1f'%(yyyy,mm,int(dd),h,min,sec)
    ret = ret.replace(' ','0')
    ret = ret.replace('A',' ')
    return ret
#print jd2gd(sys.argv[1])
