# This script calls the extinction procedure for different distances in a loop.
import numpy as np
import math 
# Define the galactic longitude and latitude
#glong = 63.0   # Example galactic longitude in degrees
#glat = -6.00   # Example galactic latitude in degrees

#glong = 12.00
#glat = -63.00

# glong = 4.5
# glat = -72.5

#glong = 67.00
#glat = -13.00

#glong = 270.0
#glat = 4.0

#glong = 270.0
#glat = -73.0
def list_ebv(glong,glat):
    
    def extin(glong,glat,dist):
        
        r0 = 8.0                                                # adopted distance of the Galactic center
        conv = np.pi/180.
        step = 0.05                                             #steps of the gas density integration to obtain column 
                                                                # density, in pc
    
        nstep = int(dist/step)
    
        if nstep == 0:
            nstep = 1
    
        #computes  trigonometric functions only once
        yproj = np.cos(glong*conv)
        xproj = np.sin(glong*conv)
        bproj = np.sin(glat*conv)
        dproj = np.cos(glat*conv)
        av = 0.0                                                #for the integration of the colunar density
    
        # declaring and puting values in the variables. The arrays will contain the
        # value of quantities like galactic radius or gas density for each
        # step along the line-of sight. If you work with other language you
        # should probably define these quantities in a loop
    
        ipas = np.arange(1,nstep+1,1)
        nel = len(ipas)
    
        dis = np.zeros(nel)
        agas = np.zeros(nel)
        x = np.zeros(nel)
        y = np.zeros(nel)
        yy = np.zeros(nel)
        r = np.zeros(nel)
        z = np.zeros(nel)
        zCO = np.zeros(nel)
        zH = np.zeros(nel)
        zc = np.zeros(nel)
        ah1 = np.zeros(nel)
        aco = np.zeros(nel)
        zmet = np.zeros(nel)
        for i in range(nel):
            dis[i] = ipas[i]*step - step    
            x[i] = (dis[i]*xproj)*dproj
            y[i] = dis[i]*yproj*dproj
            yy[i] = r0-y[i]
            r[i] = math.sqrt(x[i]*x[i]+yy[i]*yy[i])
            z[i] = dis[i]*bproj
    
            zCO[i]=0.036*math.exp(0.08*r[i])     # H2 scale-height
            zH[i] = zCO[i]*1.8                   # H1 scale-height (Guilbert 1978)
            zc[i] = 0.02                         # shift takes in to account that the sun is not 
                                                 # precisely in the galactic plane
    
            ah1[i]=0.7*math.exp(-r[i]/7.0-((1.9/r[i])**2))                                   # function that calculates the HI density
            aco[i] = 58.*math.exp(-r[i]/1.20-((3.50/r[i])**2)) + 240.*math.exp(-(r[i]**2/0.095))   # H2 density; last term is for galactic center
    
            if r[i] < 1.2:
                zmet[i] = 9.6
            
            if r[i] > 1.2 and r[i] < 9.0:
                zmet[i] = (r0/r[i])**0.5
        
            if r[i] > 9.0:
                zmet[i] = (r0/r[i])**0.1
            # this defines the metallicity correction, see section 3 of the paper
        
        ah1[0] = 0.0
        aco[0] = 0.0
        gam1=1.0
        gam2=2.0
    
        #See the final tuning (section 4.1) correction factor for interval l=120-200
    
        tune=1.
        if glong > 120 and glong < 200:
            tune=2.
    
        rs = 3.05                         #ratio between total to selective extinction
        for i in range(nel):
            agas[i] = gam1*(ah1[i]*zmet[i]*math.exp(-0.5*((z[i]-zc[i])/zH[i])**2))+gam2*aco[i]*math.exp(-0.5*((z[i]-zc[i])/zCO[i])**2)
        av = np.sum(agas)*step*3.086*.57*tune    
        ebv = av/rs
        # "total" instruction gives the sum of the array elements
        # it is equivaletn to integrate along the line-of-sight. The step is in units of kpc=
        # 3.08 *10^21 cm and the conversion factor gamma= .53 10^-21 mag cm2
    
        return ebv,av
    
    # Define the range of distances (in kpc) you want to loop through
    start_dist = 0.01  # Starting distance
    end_dist = 9.92   # Ending distance
    step_dist = 0.1   # Step size
    
    # Calculate the number of steps
    nsteps = int((end_dist - start_dist) / step_dist) + 1
    
    # Create arrays to store the results
    distances = np.zeros(nsteps)
    ebv = np.zeros(nsteps)
    av = np.zeros(nsteps)
    
    # Loop through the distances and call extinspiral
    for i in range(nsteps):
        dist = start_dist + i * step_dist
        distances[i] = dist
        EBV,av = extin(glong, glat, dist)
        ebv[i] = EBV
        #print(ebv)
    
    return ebv
    # Print the results
    # print, 'EBV Spiral'
    # for i in range(nsteps):
    #     print(ebv[i],',')
        