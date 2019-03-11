# -*- coding: utf-8 -*-
"""
Alston co-ords
53° 48' 05.14" N 53.81
02° 35' 25.53" W -02.56
"""
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import coordinates as coord
from astropy.time import Time
from astropy import units as u


observingNightMidnightUTC = Time('2019-04-25 00:00:00')
alston=EarthLocation(lat=53.81 *u.deg, lon=-02.56*u.deg, height=50*u.m)

#wasp12b values
wasp12b =SkyCoord('06 30 33 +29 40 20', unit=(u.hourangle, u.deg))
ephemeralPeriod = 1.0914203 # days
linearEphemeris=2456176.668258
transitFWHM = 24 * 0.12483 #hours


print(observingNightMidnightUTC)
epochNumber= int(round((observingNightMidnightUTC.jd - linearEphemeris)/ephemeralPeriod))
print('epoch# ='+str(epochNumber))

obsDayTimeofEphemerisMidPoint = Time(linearEphemeris + epochNumber *  ephemeralPeriod, format='jd', scale='utc')
#obsDayTimeofEphemerisEndPoint = Time(linearEphemeris + epochNumber * 1.0914203, format='jd', scale='utc')+(transitFWHM/2)*u.day
obsDayTimeofEphemerisMidPoint.format = 'iso'
print('obsDayTimeofEphemerisMidPoint=')
print(obsDayTimeofEphemerisMidPoint)

offsetFromMidnight =  ((obsDayTimeofEphemerisMidPoint - observingNightMidnightUTC)).sec/3600 #this is to align the plots



deltaFromObservingNightMidnightUTC= np.linspace(-6, 6, 100)*u.hour #to create x.axis values for the 6 hours either side of midnight
lttObservingNight =Time.light_travel_time(observingNightMidnightUTC+deltaFromObservingNightMidnightUTC, wasp12b, location=alston)
bjdhoursObservingNight =lttObservingNight.to(u.hour)
deltaFromEphemerisMidPoint= np.linspace(-(transitFWHM/2),(transitFWHM/2), 20)*u.hour #to create x.axis values for half the FHWM period either side of the transit midpoint
lttEphemeris =Time.light_travel_time(obsDayTimeofEphemerisMidPoint+deltaFromEphemerisMidPoint, wasp12b, location=alston)
bjdhoursEphemeris =lttEphemeris.to(u.hour)



print('transitFWHM/2= ' + str(transitFWHM/2))

frameOfObservingNight = AltAz(obstime=observingNightMidnightUTC+deltaFromObservingNightMidnightUTC,location=alston)
frameOfEphemeris = AltAz(obstime=obsDayTimeofEphemerisMidPoint+deltaFromEphemerisMidPoint,location=alston)

wasp12bAirmasssOfObservingNight = wasp12b.transform_to(frameOfObservingNight).secz
wasp12bAirmasssOfEphemeris = wasp12b.transform_to(frameOfEphemeris).secz



print('>>>>>Start and stop of Tranist UTC >')

print (obsDayTimeofEphemerisMidPoint+deltaFromEphemerisMidPoint[0]- bjdhoursEphemeris[0])
print (obsDayTimeofEphemerisMidPoint+deltaFromEphemerisMidPoint[19]- bjdhoursEphemeris[19])

plt.subplot(1,1,1)
plt.plot(deltaFromObservingNightMidnightUTC -bjdhoursObservingNight, wasp12bAirmasssOfObservingNight,color='gray')
plt.plot(deltaFromEphemerisMidPoint+offsetFromMidnight*u.h - bjdhoursEphemeris , wasp12bAirmasssOfEphemeris,'g',linewidth=4.0) #convert BJD to UTC
plt.xlim(-6, 6)
plt.ylim(1, 4)
plt.xlabel('Hours from '+ observingNightMidnightUTC.datetime.strftime("%d %B %Y  %I:%M%p") + ' UTC')
plt.ylabel('WASP-12 b Airmass' )

plt.show()
