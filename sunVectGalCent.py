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
theta_d = 0 # Composition dipole pointed in some direction

for dayNum in days:

	timeStamp=dayNum*24.0*3600.0+newYearTime

	dateTime = datetime.fromtimestamp(timeStamp)
	print(dayNum)

	galLocGal=SkyCoord(l=0*u.degree, b=0*u.degree, frame='galactic')

	galLocLocal=galLocGal.transform_to(AltAz(obstime=dateTime,location=Seattle))

	galDec=galLocLocal.alt.deg
	galAsc=galLocLocal.az.deg

	galTheta=90.0-galDec-labLat
	galPhi=galAsc-labLong

	#inPhase[index]=galDec
	#outPhase[index]=galAsc
	inPhase[index]=-np.cos(galDec*np.pi/180)*np.sin((galAsc+theta_d)*np.pi/180)
	outPhase[index]= -np.cos(galDec*np.pi/180)*np.cos((galAsc+theta_d)*np.pi/180)
	index+=1

output=np.column_stack((days,inPhase,outPhase))
np.savetxt('galVectMin.out',output)

plt.plot(days,inPhase)
plt.plot(days,outPhase)
plt.show()
