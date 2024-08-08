from astropy.time import Time
from astropy.coordinates import get_sun
from astropy.coordinates import SkyCoord
from astropy.coordinates import GCRS
from astropy.coordinates import CartesianRepresentation
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz
from datetime import datetime
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np

days=np.arange(1.0,366.0,60.0/(24.0*3600.0))
index=0
labLat=47.659970
labLong=-122.303400
newYearTime=1704096000
Seattle = EarthLocation(lat=labLat*u.deg,
                            lon=labLong*u.deg, height=56*u.m)
theta_d = 0 # Composition dipole pointed in some direction

## X-Vector

inPhaseX=np.zeros(len(days))
outPhaseX=np.zeros(len(days))

for dayNum in days:

	timeStamp=dayNum*24.0*3600.0+newYearTime

	dateTime = datetime.fromtimestamp(timeStamp)
	print(dayNum)

	xGeo=SkyCoord(ra=0*u.degree, dec=0*u.degree, frame='gcrs')
	xLocal=xGeo.transform_to(AltAz(obstime=dateTime,location=Seattle))

	xAlt=xLocal.alt.deg
	xAz=xLocal.az.deg

	inPhaseX[index]=-np.cos(xAlt*np.pi/180)*np.sin((xAz+theta_d)*np.pi/180)
	outPhaseX[index]= -np.cos(xAlt*np.pi/180)*np.cos((xAz+theta_d)*np.pi/180)
	index+=1

output=np.column_stack((days,inPhaseX,outPhaseX))
np.savetxt('xVectMin.out',output)

## Y-Vector

inPhaseY=np.zeros(len(days))
outPhaseY=np.zeros(len(days))

for dayNum in days:

	timeStamp=dayNum*24.0*3600.0+newYearTime

	dateTime = datetime.fromtimestamp(timeStamp)
	print(dayNum)

	yGeo=SkyCoord(ra=90*u.degree, dec=0*u.degree, frame='gcrs')
	yLocal=yGeo.transform_to(AltAz(obstime=dateTime,location=Seattle))

	yAlt=yLocal.alt.deg
	yAz=yLocal.az.deg

	inPhaseY[index]=-np.cos(yAlt*np.pi/180)*np.sin((yAz+theta_d)*np.pi/180)
	outPhaseY[index]= -np.cos(yAlt*np.pi/180)*np.cos((yAz+theta_d)*np.pi/180)
	index+=1

output=np.column_stack((days,inPhaseY,outPhaseY))
np.savetxt('yVectMin.out',output)

## Z-Vector

inPhaseZ=np.zeros(len(days))
outPhaseZ=np.zeros(len(days))

for dayNum in days:

	timeStamp=dayNum*24.0*3600.0+newYearTime

	dateTime = datetime.fromtimestamp(timeStamp)
	print(dayNum)

	zGeo=SkyCoord(ra=0*u.degree, dec=90*u.degree, frame='gcrs')
	zLocal=zGeo.transform_to(AltAz(obstime=dateTime,location=Seattle))

	zAlt=zLocal.alt.deg
	zAz=zLocal.az.deg

	inPhaseZ[index]=-np.cos(zAlt*np.pi/180)*np.sin((zAz+theta_d)*np.pi/180)
	outPhaseZ[index]= -np.cos(zAlt*np.pi/180)*np.cos((zAz+theta_d)*np.pi/180)
	index+=1

output=np.column_stack((days,inPhaseZ,outPhaseZ))
np.savetxt('zVectMin.out',output)
