pro extinspiral, glong, glat, dist, ebvspiral, avspiral

; EXTINSPIRAL.PRO tested on feb 1, 2005
; This program computes extinction between Sun and a point with given distance
; and galactic coordinates
; it requires functions SPIRAL2 and ARMDENS in the same folder to work

FORWARD_FUNCTION spiral2, armdens

r0=8.27          ; adopted distance of Sun to galactic center, in kpc
conv=!pi/180.
step  = 0.05     ; steps in kpc for integration of density to obtain column 
                 ; density, in pc

nstep=fix(dist/step)  ; number of elements of arrays that describe extinction

if nstep eq 0 then nstep=1

;computes  trigonometric functions only once for a given line-of-sight

yproj=cos(glong*conv)
xproj=sin(glong*conv)
bproj=sin(glat*conv)
dproj=cos(glat*conv)
;av=0.0                  ;for the integration of the colunar density
along=glong*conv          ;longitude in radians


;declaring arrays and puting initial values. The arrays will contain the
;value of quantities like galactic radius or gas density for each step along the line-of sight
;if you work with other language you should probably define these quantities in a loop
 dis= make_array(nstep,/float,value=0.0)
x  = make_array(nstep,/float,value=0.0)
y  = make_array(nstep,/float,value=0.0)
yy = make_array(nstep,/float,value=0.0)
r  = make_array(nstep,/float,value=0.0)
z  = make_array(nstep,/float,value=0.0)
zCO= make_array(nstep,/float,value=0.0)
zH = make_array(nstep,/float,value=0.0)
ah1= make_array(nstep,/float,value=0.0)
aco= make_array(nstep,/float,value=0.0)
zmet=make_array(nstep,/float,value=0.0)
agas=make_array(nstep,/float,value=0.0)
escz=make_array(nstep,/float,value=0.0)
esczH1=make_array(nstep,/float,value=0.0)
densh1=make_array(nstep,/float,value=0.0)
densh2=make_array(nstep,/float,value=0.0)
h2bckg=make_array(nstep,/float,value=0.0)
h1bckg=make_array(nstep,/float,value=0.0)
toth2=make_array(nstep,/float,value=0.0)
spidens=make_array(nstep,/float,value=0.0)
acol= make_array(nstep,/float,value=0.0)
av = make_array(nstep,/float,value=0.0)
ebv = make_array(nstep,/float,value=0.0)
axisym = make_array(nstep,/float,value=0.0)
ebvaxysim = make_array(nstep,/float,value=0.0)
ipas=findGen(nstep)/1 +1  ; generates an array with a sequence of numbers, used as index for
                          ; distance along line-of-sight
NUMEL=n_elements(ipas)

dis=ipas*step  - step           ;distance along l-o-s
x=(dis*xproj)*dproj
y=dis*yproj*dproj
yy=r0-y
r=sqrt(x*x+yy*yy)         ;galactic radius along l-o-s
z=dis*bproj

zCO=0.036*exp(0.08*r)     ;H2 scale-height
zH = zco*1.8              ;H1 scale-height (Guilbert 1978)
zc = 0.02                 ;shift takes in to account that the sun is not precisely in the galactic plane

ah1=0.7*exp(-r/7.0-((1.9/r)^2))  ;function that calculates the HI density
aco = 58.*exp(-r/1.20-((3.5/r)^2)) + 240.*exp(-(r^2/0.095)) ;H2 density; last term is for galactic center region

ah1[0] = 0.0
aco[0] = 0.0

for i=0,NUMEL-1 do begin
    if r[i] le 1.2 then  zmet[i] = 9.6
    if r[i] gt 1.2 and r[i] le 9.0 then zmet[i] = (r0/r[i])^0.5
    if r[i] gt 9.0 then  zmet[i] = (r0/r[i])^0.1
endfor
         ; this defines the metallicity correction, see section 3 of the paper


escz=exp(-0.50*((z-zc)/zCO)^2)     ;density variation due to varittion in distance xz to the plae
esczHI=exp(-0.50*((z-zc)/zH)^2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;next we create the arrays of points representing an arm, with function spiral2
;then we calculate the contribution of the arm to the column density using function armdens
; there is a list of arms for HI gas distribution and another for H2 (see paper)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    h2bckg= escz*aco  ;background density along l-o-s, corresponding to axissimetric model


    b1h2=spiral2(2.70,-20.0,7.0,0.0,380);                ;this creates an arm, see list in the paper
    spidens=armdens(b1h2,along,0.1,0.125,nstep)*h2bckg   ; creates an array containing densities due to that arm
    densh2=densh2+spidens                                ;contribution of each arm is added


    bb4h2=spiral2(2.4,165.,12.3,0.00,360);  1
    spidens=armdens(bb4h2,along,3.0,0.135,nstep)*h2bckg
    densh2=densh2+spidens

    b41h2=spiral2(6.5,353.,13.8,0.1,120);*  3
    spidens=armdens(b41h2,along,6.,0.08,nstep)*h2bckg
    densh2=densh2+spidens

    b42h2=spiral2(3.65, 30.,13.5,0.05,150);4
    spidens=armdens(b42h2,along,10.,0.140,nstep)*h2bckg
    densh2=densh2+spidens

    b43h2=spiral2(2.50, 20.,7.5,0.00,345);5
    spidens=armdens(b43h2,along,0.5,0.06,nstep)*h2bckg
    densh2=densh2+spidens

    b44h2=spiral2(2.9, -75.,11.7,0.05,180);6
    spidens=armdens(b44h2,along,1.0,0.080,nstep)*h2bckg
    densh2=densh2+spidens

    b55h2=spiral2(6.8,-30.0,6.0,0.0,80);7 *
    spidens=armdens(b55h2,along,24.,0.040,nstep)*h2bckg
    densh2=densh2+spidens

    b88h2=spiral2(10.10,-120.0,7.10,0.0,140);8 *
    spidens=armdens(b88h2,along,0.80,0.125,nstep)*h2bckg
    densh2=densh2+spidens

    b99h2=spiral2(8.40,-140.0,6.95,0.0,140);9 *
    spidens=armdens(b99h2,along,2.0,0.125,nstep)*h2bckg
    densh2=densh2+spidens

    bb3h2=spiral2(7.5,-30.0,8.0,0.0,40.); 10
    spidens=armdens(bb3h2,along,4.,0.08,nstep)*h2bckg
    densh2=densh2+spidens

    badh2=spiral2(7.9,-70.,10.0,0.0,130);11
    spidens=armdens(badh2,along,10.0,0.145,nstep)*h2bckg
    densh2=densh2+spidens

    b77h2=spiral2(5.70, 5.,10.0,-0.15,150);*
    spidens=armdens(b77h2,along,8.0,0.135,nstep)*h2bckg
    densh2=densh2+spidens

    b71h2=spiral2(5.1,345.,55.,0.0,15);* 13
    spidens=armdens(b71h2,along,10.0,0.125,nstep)*h2bckg
    densh2=densh2+spidens

    b81h2=spiral2(6.20, 185.,6.5,0.0,110);*  14
    spidens=armdens(b81h2,along,30.0,0.125,nstep)*h2bckg

    densh2=densh2+spidens


toth2=1.0*densh2+0.1*h2bckg  ;contribution of H2 arms plus a small amount of background for interarms


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;HI arms

h1bckg= esczHI*zmet*ah1         ; the axis-symmetric contribution due to HI gas along the line-o-s
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    b1hi =spiral2(2.35,-31,6.5,0.0,240)
    spidens=armdens(b1hi,along,0.1,0.205,nstep)*h1bckg
    densh1=densh1+spidens

    b91hi=spiral2(3.80,195.,6.70,0.0,330)
    spidens=armdens(b91hi,along,0.3,0.13,nstep)*h1bckg
    densh1=densh1+spidens

    b2hi=spiral2(2.57,155.,7.1,0.0,200)
    spidens=armdens(b2hi,along,0.02,0.15,nstep)*h1bckg
    densh1=densh1+spidens

    b81hi=spiral2(3.8,-14.,7.5,0.0,245)
    spidens=armdens(b81hi,along,0.2,0.175,nstep)*h1bckg
    densh1=densh1+spidens

    b41hi=spiral2(7.50,-2.5,6.35,0.1,182)
    spidens=armdens(b41hi,along,0.3,0.125, nstep)*h1bckg
    densh1=densh1+spidens

    b42hi=spiral2(7.65, 163.,7.5,0.1,159)        ;
    spidens=armdens(b42hi,along,1.0,0.180,nstep)*h1bckg
    densh1=densh1+spidens

    b43hi=spiral2(5.9,222.,5.3,0.21,279)
    spidens=armdens(b43hi,along,0.3,0.25,nstep)*h1bckg
    densh1=densh1+spidens

    b44hi=spiral2(9.05, 345.,-22.,-1.8,24)  ;testado
    spidens=armdens(b44hi,along,4.0,0.135,nstep)*h1bckg
    densh1=densh1+spidens

    bb43hi=spiral2(6.44, 225.,11.55,-0.10,110)
    spidens=armdens(bb43hi,along,0.2,0.22,nstep)*h1bckg
    densh1=densh1+spidens

    bb44hi=spiral2(7.36,-8.4,7.75,0.0,199)
    spidens=armdens(bb44hi,along,2.5,0.18,nstep)*h1bckg
    densh1=densh1+spidens

    b51hi=spiral2(3.95,190.,9.4,0.0,170)
    spidens=armdens(b51hi,along,0.2,0.175,nstep)*h1bckg
    densh1=densh1+spidens

    b201hi=spiral2(8.60,10.50,45.40,0.0,69)
    spidens=armdens(b201hi,along,0.5,0.210,nstep)*h1bckg
    densh1=densh1+spidens

    bb7hi=spiral2(7.3,330.0,12.0,0.0,031)
    spidens=armdens(bb7hi,along,8.,0.25,nstep)*h1bckg
    densh1=densh1+spidens

    b71hi=spiral2(11.4,315.,-40.,0.0,31)
    spidens=armdens(b71hi,along,2.0,0.135,nstep)*h1bckg
    densh1=densh1+spidens

    bb13hi=spiral2(10.1,-10,11.1,-0.3,48) ;viol
    spidens=armdens(bb13hi,along,2.0,0.190,nstep)*h1bckg
    densh1=densh1+spidens

    bb211hi=spiral2(5.0,-45.0,9.0,-0.3,25) ;viol
    spidens=armdens(bb211hi,along,3.0,0.125,nstep)*h1bckg
    densh1=densh1+spidens


    bb2hi=spiral2(11.1,-38.,11.1,-0.05,65)
    spidens=armdens(bb2hi,along,1.0,0.190,nstep)*h1bckg
    densh1=densh1+spidens

    b101hi=spiral2(8.6,25.4,10.5,0.0,89)
    spidens=armdens(b101hi,along,1.0,0.285,nstep)*h1bckg
    densh1=densh1+spidens



densh1= densh1 + 0.75*h1bckg  ;add the small contribution of HI arms to large axis-symmetric background

gam1=1.0
gam2=2.0

;See the final tuning (section 4.) for correction between l=120 and l=200

finetune=1.
if (glong ge 119.) and (glong le 200.) then finetune=2.

     dentot=gam1*densh1+gam2*toth2  ;this is the extinction distribution along l-o-s  with spiral arms

    axisym=(gam1*h1bckg + gam2*h2bckg)*finetune  ;this is axysymmetric model if you want it

    acol=total(dentot, /cumulative) ; produce array containing integrated column density
    ; the cumulative instruction replace each element of an array with the sum of the elements up to that point

   ;conversion factor below includes dis(kpc)=3.086e21 cm and av(mag cm2)= 5.7e-22 NH see paper

    fatconv=step*3.086*.57  ; step(kpc)*3.086e21/5.7 e-22
    av=acol*fatconv
    rs = 3.05
    ebv = av/rs
;up to this moment you have arrays containing extinction as a function of distance along the l-o-s
;to print the extinction to the distance you entered, just take the last element

avspiral=av[numel-1]
ebvspiral=avspiral/rs

floating_point_underflow = 32
status = Check_Math()         ; Get status and reset accumulated math error register.
IF(status AND NOT floating_point_underflow) NE 0 THEN $
    Message, 'IDL Check_Math() error: ' + StrTrim(status, 2)


end


