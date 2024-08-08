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

#days=np.arange(200.0,311.0,1.0/(24.0*3600.0))
days=np.arange(1.0,366.0,60.0/(24.0*3600.0))
inPhase=np.zeros(len(days))
outPhase=np.zeros(len(days))
index=0
labLat=47.659970
labLong=-122.303400
newYearTime=1704096000
Seattle = EarthLocation(lat=labLat*u.deg,
                            lon=labLong*u.deg, height=56*u.m)
theta_d = -31 # Composition dipole pointed in some direction

for dayNum in days:

	timeStamp=dayNum*24.0*3600.0+newYearTime

	dateTime = datetime.fromtimestamp(timeStamp)
	print(dayNum)

	sunLocGeo=get_sun(Time(dateTime))

	sunLocLocal=sunLocGeo.transform_to(AltAz(obstime=dateTime,location=Seattle))

	sunDec=sunLocLocal.alt.deg
	sunAsc=sunLocLocal.az.deg

	sunTheta=90.0-sunDec-labLat
	sunPhi=sunAsc-labLong

	#inPhase[index]=sunDec
	#outPhase[index]=sunAsc
	inPhase[index]=-np.cos(sunDec*np.pi/180)*np.sin((sunAsc+theta_d)*np.pi/180)
	outPhase[index]= -np.cos(sunDec*np.pi/180)*np.cos((sunAsc+theta_d)*np.pi/180)
	index+=1

output=np.column_stack((days,inPhase,outPhase))
np.savetxt('sunVectMin85-97.out',output)

plt.plot(days,inPhase)
plt.plot(days,outPhase)
plt.show()
