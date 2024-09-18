# EXTINSPIRAL.PRO tested on feb 1, 2005
# This program computes extinction between Sun and a point with given distance
# and galactic coordinates
# it requires functions SPIRAL2 and ARMDENS in the same folder to work
import numpy as np
import math

# glong = 63.0     #118
# glat = -6.0     #67.9

def list_ebv(glong,glat):        
    
    def extin_spiral(glong,glat,dist):
        
        r0=8.00          # adopted distance of Sun to galactic center, in kpc  --> r0=8.27
        conv=np.pi/180.
        step  = 0.05     # steps in kpc for integration of density to obtain column 
                         # density, in pc
    
        nstep=int(dist/step)  # number of elements of arrays that describe extinction
        
        
        if nstep == 0:
            nstep=1
        
        #computes  trigonometric functions only once for a given line-of-sight
        
        yproj=math.cos(glong*conv)
        xproj=math.sin(glong*conv)
        bproj=math.sin(glat*conv)
        dproj=math.cos(glat*conv)
        #av=0.0                  #for the integration of the colunar density
        along=glong*conv          #longitude in radians
        
        ipas = np.arange(1,nstep+1,1)
        NUMEL = len(ipas)
        
        # def armdens(arm,glong,denmax,sig,nstep):
        
        #     #This function creates an array containing "intensities" along a line-of-sight
        #     #produced by the presence of a spiral arm. The number of elements of this array
        #     #is given by NSTEP , defined in the main program (distance divided by step-0.05 kpc)
        #     #ARM is the array describing the arm, created by the function spiral.
        #     #The array contains x,y, radius and polar angle theta of all the points of the arm
        #     #(refered to the Galactic center)
        #     #long is the longitude (in radians) of the line-of-sight. The arm is supposed to
        #     #have a Gaussian "density" distribution, with peak value DENMAX and width SIG
        #     #The first step is to verify if there is an intersection of the arm with the line-of-sight
        #     #  ;all the calculations refer to the gas density in the Galactic plane; in the main program
        #     #  ;that calls this routine the effect of Galactic latitude different from zero is taken into
        #     #  ;account
        
        #     #*********************************************************************
        #     #next find the number of elements contained in array arm
        #     arm = arm
        #     nel = int(n_elements(arm)/4)
        
        #     step = 0.05
        #     #dist = step*findgen(nstep) #creates an array with distances along the line-o-s, in steps of 0.05 kpc
        #     dist = np.arange(0,nstep,0.05)
            
        #     # reserves memory for arrays
        #     rgal = np.zeros(nstep)
        #     dr = np.zeros(nstep)
        #     dens = np.zeros(nstep)
        
        #     r0 = 7.5   #adopted distance of Galactic center, in kpc
        #     conv = np.pi/180
        #     for i in range(nstep):
        #         rgal[i] = math.sqrt(dist[i]**2 + r0**2 -2*dist[i]*r0*math.cos(glong))
        #         #creates array with  galactic radius of each point along the line-of-sight
        
        #     #****** We next compute the Galactic longitude of each point of the arm
        #     armgl = np.zeros(nel)
        #     for i in range(1,nel+1,1):
        #         if abs(arm[i-1,1] -r0) < 0.0001 and arm[i-1,0] > 0.:
        #             armgl[i-1] = np.pi/2.   #this mean y very close to zero
        #         if abs(arm[i-1,1] -r0) < 0.0001 and arm[i-1,0] < 0.:
        #             armgl[i-1] = -np.pi/2
        #         if (arm[i-1,1] -r0) < -0.0001:
        #             angl = math.tan(arm[i-1,0]/abs(r0-arm[i-1,1]))    #use the x and y coordinates with respect to Sun to find longitude
        #             armgl[i-1] = angl
                
        #         if ((arm[i-1,1] -r0) > 0.0001):
        #             angl = math.tan(arm[i-1,0]/abs(arm[i-1,1]-r0)) #tenho duvidas
        #             armgl[i-1] = np.pi-angl
                
        #         if armgl[i-1] < 0.:
        #             armgl[i-1] = 2*np.pi+ armgl[i-1]
        
        #     #****** We next look for the point where the arm crosses the lin-of-sight
        #     ncross = np.zeros(4)
        #     dlong = np.zeros(nel)
        #     dlong = (armgl - glong)    ;difference in longitude between points of arm and l-o-sight
        
        #     jj = 0
        #     for i in range(0,nel-2,1):
                
        #         if dlong[i] == 0.:
        #             ncross[jj] = i
        #         if dlong[i]*dlong[i+1] < 0. and abs(armgl[i]-armgl[i+1]) < 0.5:
        #             ncross[jj] = i  #this deals with fact that a l-o-s can cross an arm several times (some armas make several turns)
        #             jj+=1
        #        #if dlong changes  sign, the arm is crossing the line-of-sight.
        #        # we keep in memory ncross, the number of the element of the arm  array where arm is crossed by l-o-s
        
        #        #ncr refers to the points along the arms; kbrac to the points along the line-of-sight
        #     ncr = ncross[0]
        #     if ncr > 1:
        #         dsun=sqrt(arm[ncr,0]^2 + (arm[ncr,1]-r0)^2)
        #     karm=fix(dsun/step)
        #      ;karm is the index of the point of the l-o-s where it crosses the arm, will be center of Gaussian distribution
        
        #        ;** calculation of angle between arm and l-o-s, in order to obtain the width of the Gaussian
        #        ;density distribution
        #           alph=atan((arm[ncr+1,1]-arm[ncr-1,1])/(arm[ncr+1,0]-arm[ncr-1,0]))
        #        ;delta Y/delta X gives the angle of arm with respect to horizontal axis
        #           beta=glong-!pi/2
        #           factor=1./abs(sin(glong+alph))
        #         if factor ge 10 then factor=10.
        #         ;this limits the width of the Gaussian to a factor 10 larger than when l-o-s crosses perpendicular
        
        #         sig2=2*(factor*sig/.05)^2
        
                 #armgl=armgl/conv
        def armdens(arm, glong, denmax, sig, nstep):
            """
            This function creates an array containing "intensities" along a line-of-sight
            produced by the presence of a spiral arm. The number of elements of this array
            is given by nstep.
            
            Parameters:
            arm    : numpy array containing x, y, radius, and polar angle theta of all the points of the arm
            glong  : Galactic longitude (in radians) of the line-of-sight
            denmax : Peak value of the Gaussian "density" distribution
            sig    : Width of the Gaussian "density" distribution
            nstep  : Number of steps along the line-of-sight
        
            Returns:
            dens   : Array with density along the line-of-sight
            """
            arm = np.array(arm)
            nel = int(len(arm) / 4)
        
            step = 0.05
            dist = step + step * np.arange(nstep)  # Creates an array with distances along the line-of-sight, in steps of 0.05 kpc
            # print(dist)
            # exit()
            rgal = np.zeros(nstep, dtype=float)
            dr = np.zeros(nstep, dtype=float)
            dens = np.zeros(nstep, dtype=float)
        
            r0 = 8.0     #7.5  # Adopted distance of Galactic center, in kpc
            conv = np.pi / 180
            rgal = np.sqrt(dist**2 + r0**2 - 2 * dist * r0 * np.cos(glong))  # Galactic radius of each point along the line-of-sight
        
            armgl = np.zeros(nel, dtype=float)
        
            for i in range(1, nel + 1):
                if abs(arm[i - 1, 1] - r0) < 0.0001:
                    if arm[i - 1, 0] > 0:
                        armgl[i - 1] = np.pi / 2
                    else:
                        armgl[i - 1] = -np.pi / 2
                if arm[i - 1, 1] - r0 < -0.0001:
                    angl = np.arctan(arm[i - 1, 0] / abs(r0 - arm[i - 1, 1]))
                    armgl[i - 1] = angl
                if arm[i - 1, 1] - r0 > 0.0001:
                    angl = np.arctan(arm[i - 1, 0] / abs(arm[i - 1, 1] - r0))
                    armgl[i - 1] = np.pi - angl
                
                if armgl[i - 1] < 0:
                    armgl[i - 1] = 2 * np.pi + armgl[i - 1]
        
            ncross = np.zeros(4, dtype=int)
            dlong = armgl - glong
        
            jj = 0
            for i in range(nel - 2):
                if dlong[i] == 0:
                    ncross[jj] = i
                if dlong[i] * dlong[i + 1] < 0 and abs(armgl[i] - armgl[i + 1]) < 0.5:
                    ncross[jj] = i
                    jj = jj+ 1
        
            # First crossing
            ncr = ncross[0]
            if ncr > 1:
                dsun = np.sqrt(arm[ncr, 0]**2 + (arm[ncr, 1] - r0)**2)
                karm = int(dsun / step)
        
                alph = np.arctan((arm[ncr + 1, 1] - arm[ncr - 1, 1]) / (arm[ncr + 1, 0] - arm[ncr - 1, 0]))
                beta = glong - np.pi / 2
                factor = 1.0 / abs(np.sin(glong + alph))
                factor = min(factor, 10)
        
                sig2 = 2 * (factor * sig / 0.05)**2
        
                for k in range(nstep):
                    karmdelt = abs(k - karm)
                    if karmdelt > 15:
                        dens[k] = 0.0
                    else:
                        dens[k] = denmax * np.exp(-((k - karm)**2) / sig2)
        
            # Second crossing (if any)
            ncr = ncross[1]
            if ncr > 1:
                dsun = np.sqrt(arm[ncr, 0]**2 + (arm[ncr, 1] - r0)**2)
                karm = int(dsun / step)
        
                alph = np.arctan((arm[ncr + 1, 1] - arm[ncr - 1, 1]) / (arm[ncr + 1, 0] - arm[ncr - 1, 0]))
                beta = glong - np.pi / 2
                gam = abs(glong + alph)
                if gam >= 0.1:
                    factor = 1.0 / np.sin(gam)
                else:
                    factor = 10.0
        
                sig2 = 2 * (factor * sig / 0.05)**2
        
                for k in range(nstep):
                    karmdelt = abs(k - karm)
                    if karmdelt > 15:
                        dens[k] = 0.0
                    else:
                        dens[k] = denmax * np.exp(-((k - karm)**2) / sig2)
        
            return dens
        
        def spiral2(r0,fase,pitch0,dpitch,npont):
            ''' this function creates a 2-dimensional array containing NPONT lines, one for each point of the arm
            ; the columns are X,Y,R,THETA coordinates of each point (all referred to the Galactic center)
            ; parameters are: r0= initial radius of the arm , fase= initial azimuthal angle (direction of Sun is fase=0)
            ;pitch0= initial pitch angle of the arms (angle with respect to circles). for most arms this is constant,
            ; but we introduced the possibility of variable (steadily increasing or decreasing) pitch angle.
            ;dpitch is the rate of variation
            '''
            pitch = np.zeros(npont)
            r = np.zeros(npont)
            teta = np.zeros(npont)
            x = np.zeros(npont)
            y = np.zeros(npont)
        
            conv = np.pi/180
            cte = math.tan(pitch0*conv)
            tet = np.arange(npont)*conv  # Creates an array with theta for each point of the arm
            faserad = fase*conv
            teta = tet+faserad
            r[0] = r0
            pitch[0] = pitch0
        
            for j in range(1,npont,1):
                pitch[j]=pitch0 + j*dpitch      #the pitch angle can vary!
                r[j] = r[j-1] *(1.+ math.tan(pitch[j-1]*conv)*conv)
            ''';above, for each radius, the next radius is the old one plus the increase in radius,
            ;which is equal to the delta angle (one degree=conv)times the tangent of the pitch angle
            ;r=r0*EXP(teta*cte) is the usual equation of an arm; we do not use it to allow variable pitch
            '''
            
            for i in range (npont):
                x[i] = -r[i]*math.sin(teta[i])
                y[i] = r[i]*math.cos(teta[i])
                
            spi = np.vstack([x, y, r, teta]).T
        
            return spi   #returns the arm array to the main program
        
        dis = np.zeros(nstep,dtype = float)
        agas = np.zeros(nstep,dtype = float)
        x = np.zeros(nstep,dtype = float)
        y = np.zeros(nstep,dtype = float)
        yy = np.zeros(nstep,dtype = float)
        r = np.zeros(nstep,dtype = float)
        z = np.zeros(nstep,dtype = float)
        zCO = np.zeros(nstep,dtype = float)
        zH = np.zeros(nstep,dtype = float)
        zc = np.zeros(nstep,dtype = float)
        ah1 = np.zeros(nstep,dtype = float)
        aco = np.zeros(nstep,dtype = float)
        zmet = np.zeros(nstep,dtype = float)
        escz = np.zeros(nstep,dtype = float)
        esczH1 = np.zeros(nstep,dtype = float)
        densh1 = np.zeros(nstep,dtype = float)
        densh2 = np.zeros(nstep,dtype = float)
        h2bckg = np.zeros(nstep,dtype = float)
        h1bckg = np.zeros(nstep,dtype = float)
        toth2 = np.zeros(nstep,dtype = float)
        spidens = np.zeros(nstep,dtype = float)
        acol = np.zeros(nstep,dtype = float)
        av = np.zeros(nstep,dtype = float)
        ebv = np.zeros(nstep,dtype = float)
        axisym = np.zeros(nstep,dtype = float)
        ebvaxysim = np.zeros(nstep,dtype = float)
        dentot = np.zeros(nstep,dtype = float)
        
        for i in range(NUMEL):
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
            escz[i] = math.exp(-0.50*((z[i]-zc[i])/zCO[i])**2)     #density variation due to varittion in distance xz to the plae
            esczH1[i] = math.exp(-0.50*((z[i]-zc[i])/zH[i])**2)
            
            if r[i] < 1.2:
                zmet[i] = 9.6
            
            if r[i] > 1.2 and r[i] < 9.0:
                zmet[i] = pow((r0/r[i]),0.5)
            
            if r[i] > 9.0:
                zmet[i] = pow((r0/r[i]),0.1)
            # this defines the metallicity correction, see section 3 of the paper
            
        ah1[0] = 0.0
        aco[0] = 0.0
        gam1 = 1.0
        gam2 = 2.0
        
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        #next we create the arrays of points representing an arm, with function spiral2
        #then we calculate the contribution of the arm to the column density using function armdens
        #there is a list of arms for HI gas distribution and another for H2 (see paper)
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        # for i in range(NUMEL):
            
        h2bckg = escz*aco  #background density along l-o-s, corresponding to axissimetric model
        
        
        b1h2 = spiral2(2.70,-20.0,7.0,0.0,380)#                #this creates an arm, see list in the paper
        spidens = armdens(b1h2,along,0.1,0.125,nstep)*h2bckg   # creates an array containing densities due to that arm
        densh2 = densh2 + spidens                                #contribution of each arm is added
        
        
        bb4h2 = spiral2(2.4,165.,12.3,0.00,360)#  1
        spidens = armdens(bb4h2,along,3.0,0.135,nstep)*h2bckg[i]
        densh2 = densh2+spidens
        
        b41h2 = spiral2(6.5,353.,13.8,0.1,120)#*  3
        spidens = armdens(b41h2,along,6.,0.08,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b42h2 = spiral2(3.65, 30.,13.5,0.05,150)#4
        spidens = armdens(b42h2,along,10.,0.140,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b43h2 = spiral2(2.50, 20.,7.5,0.00,345)#5
        spidens = armdens(b43h2,along,0.5,0.06,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b44h2 = spiral2(2.9, -75.,11.7,0.05,180)#6
        spidens = armdens(b44h2,along,1.0,0.080,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b55h2 = spiral2(6.8,-30.0,6.0,0.0,80)#7 *
        spidens = armdens(b55h2,along,24.,0.040,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b88h2 = spiral2(10.10,-120.0,7.10,0.0,140)#8 *
        spidens = armdens(b88h2,along,0.80,0.125,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b99h2 = spiral2(8.40,-140.0,6.95,0.0,140)#9 *
        spidens = armdens(b99h2,along,2.0,0.125,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        bb3h2 = spiral2(7.5,-30.0,8.0,0.0,40)# 10
        spidens = armdens(bb3h2,along,4.,0.08,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        badh2 = spiral2(7.9,-70.,10.0,0.0,130)#11
        spidens = armdens(badh2,along,10.0,0.145,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b77h2 = spiral2(5.70, 5.,10.0,-0.15,150)#*
        spidens = armdens(b77h2,along,8.0,0.135,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b71h2 = spiral2(5.1,345.,55.,0.0,15)#* 13
        spidens = armdens(b71h2,along,10.0,0.125,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
        b81h2 = spiral2(6.20, 185.,6.5,0.0,110)#*  14
        spidens = armdens(b81h2,along,30.0,0.125,nstep)*h2bckg[i]
        densh2 = densh2 + spidens
        
         
        toth2=1.0*densh2+0.1*h2bckg  #contribution of H2 arms plus a small amount of background for interarms
        
        
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        #HI arms
        
        h1bckg = esczH1*zmet*ah1         # the axis-symmetric contribution due to HI gas along the line-o-s
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        b1hi = spiral2(2.35,-31,6.5,0.0,240)
        spidens = armdens(b1hi,along,0.1,0.205,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b91hi = spiral2(3.80,195.,6.70,0.0,330)
        spidens = armdens(b91hi,along,0.3,0.13,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b2hi = spiral2(2.57,155.,7.1,0.0,200)
        spidens = armdens(b2hi,along,0.02,0.15,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b81hi = spiral2(3.8,-14.,7.5,0.0,245)
        spidens = armdens(b81hi,along,0.2,0.175,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b41hi = spiral2(7.50,-2.5,6.35,0.1,182)
        spidens = armdens(b41hi,along,0.3,0.125, nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b42hi = spiral2(7.65, 163.,7.5,0.1,159)        ;
        spidens = armdens(b42hi,along,1.0,0.180,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b43hi = spiral2(5.9,222.,5.3,0.21,279)
        spidens = armdens(b43hi,along,0.3,0.25,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b44hi = spiral2(9.05, 345.,-22.,-1.8,24)  #testado
        spidens = armdens(b44hi,along,4.0,0.135,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        bb43hi = spiral2(6.44, 225.,11.55,-0.10,110)
        spidens = armdens(bb43hi,along,0.2,0.22,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        bb44hi = spiral2(7.36,-8.4,7.75,0.0,199)
        spidens = armdens(bb44hi,along,2.5,0.18,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b51hi = spiral2(3.95,190.,9.4,0.0,170)
        spidens = armdens(b51hi,along,0.2,0.175,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b201hi = spiral2(8.60,10.50,45.40,0.0,69)
        spidens = armdens(b201hi,along,0.5,0.210,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        bb7hi = spiral2(7.3,330.0,12.0,0.0,31)
        spidens = armdens(bb7hi,along,8.,0.25,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b71hi = spiral2(11.4,315.,-40.,0.0,31)
        spidens = armdens(b71hi,along,2.0,0.135,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        bb13hi = spiral2(10.1,-10,11.1,-0.3,48) #viol
        spidens = armdens(bb13hi,along,2.0,0.190,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        bb211hi = spiral2(5.0,-45.0,9.0,-0.3,25) #viol
        spidens = armdens(bb211hi,along,3.0,0.125,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        bb2hi = spiral2(11.1,-38.,11.1,-0.05,65)
        spidens = armdens(bb2hi,along,1.0,0.190,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        b101hi = spiral2(8.6,25.4,10.5,0.0,89)
        spidens = armdens(b101hi,along,1.0,0.285,nstep)*h1bckg[i]
        densh1 = densh1 + spidens
        
        
        
        densh1 = densh1 + 0.75*h1bckg  #add the small contribution of HI arms to large axis-symmetric background
        
        #See the final tuning (section 4.) for correction between l=120 and l=200
        
        finetune = 1
        if glong > 119 and glong < 200:
            finetune = 2
        
        dentot = gam1*densh1 + gam2*toth2  #this is the extinction distribution along l-o-s  with spiral arms
        
        axisym = (gam1*h1bckg + gam2*h2bckg)*finetune  #this is axysymmetric model if you want it
        
        acol = np.sum(dentot) # produce array containing integrated column density
        # the cumulative instruction replace each element of an array with the sum of the elements up to that point
        
        #conversion factor below includes dis(kpc)=3.086e21 cm and av(mag cm2)= 5.7e-22 NH see paper
    
        fatconv=step*3.086*.57  # step(kpc)*3.086e21/5.7 e-22
        av=acol*fatconv
        rs = 3.05
        ebv = av/rs
        #up to this moment you have arrays containing extinction as a function of distance along the l-o-s
        #to print the extinction to the distance you entered, just take the last element
        
        avspiral=av
        ebvspiral=avspiral/rs
        
        return (ebvspiral)
        
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
        EBV = extin_spiral(glong, glat, dist)
        ebv[i] = EBV
        #print(ebv)
    return ebv
    
# print(extin_spiral_list(glong,glat))
    
    
    # # Print the results
    # # print, 'EBV Spiral'
    # for i in range(nsteps):
    #     print(ebv[i],',')