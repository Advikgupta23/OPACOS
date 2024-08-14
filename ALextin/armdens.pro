FUNCTION armdens, arm, glong, denmax,  sig, nstep
;This function creates an array containing "intensities" along a line-of-sight
;produced by the presence of a spiral arm. The number of elements of this array
;is given by NSTEP , defined in the main program (distance divided by step-0.05 kpc)
;ARM is the array describing the arm, created by the function spiral.
;The array contains x,y, radius and polar angle theta of all the points of the arm
;(refered to the Galactic center)
;long is the longitude (in radians) of the line-of-sight. The arm is supposed to
;have a Gaussian "density" distribution, with peak value DENMAX and width SIG
;The first step is to verify if there is an intersection of the arm with the line-of-sight
   ;all the calculations refer to the gas density in the Galactic plane; in the main program
   ;that calls this routine the effect of Galactic latitude different from zero is taken into
   ;account

;*********************************************************************
;next find the number of elements contained in array arm
arm=arm
nel=fix(n_elements(arm)/4)

  step=0.05
  dist=step*findgen(nstep) ;creates an array with distances along the line-o-s, in steps of 0.05 kpc

  ; reserves memory for arrays
  rgal=make_array(nstep, /float)
  dr=make_array(nstep, /float)
  dens=make_array(nstep, /float, value=0.)

  r0=7.5   ;adopted distance of Galactic center, in kpc
  conv=!pi/180
  rgal=sqrt(dist^2 + r0^2 -2*dist*r0*cos(glong))
   ;creates array with  galactic radius of each point along the line-of-sight

;****** We next compute the Galactic longitude of each point of the arm
armgl=make_array(nel, /float)
for i=1,nel do begin
 if abs(arm[i-1,1] -r0) lt 0.0001 and arm[i-1,0] gt 0. then armgl[i-1]=!pi/2.;this mean y very close to zero
 if abs(arm[i-1,1] -r0) lt 0.0001 and arm[i-1,0] lt 0. then armgl[i-1]=-!pi/2
 if (arm[i-1,1] -r0) lt -0.0001 then begin
    angl=atan(arm[i-1,0]/abs(r0-arm[i-1,1]))    ;use the x and y coordinates with respect to Sun to find longitude
    armgl[i-1]=angl
 endif
 if ((arm[i-1,1] -r0) gt 0.0001) then begin
   angl=atan(arm[i-1,0]/abs(arm[i-1,1]-r0)) ;tenho duvidas
              armgl[i-1]=!pi - angl
  endif
   if armgl[i-1] lt 0. then armgl[i-1]=2*!PI+ armgl[i-1]
 endfor
 ;armgl=armgl/conv

;****** We next look for the point where the arm crosses the lin-of-sight
ncross=make_array(4, /integer, value=0)
dlong=make_array(nel, /float)
dlong=(armgl - glong)    ;difference in longitude between points of arm and l-o-sight

    jj=0
   for i=0,nel-3 do begin
     if dlong[i] eq 0. then ncross[jj]=i
     if dlong[i]*dlong[i+1] lt 0. and abs(armgl[i]-armgl[i+1]) lt 0.5 then begin
        ncross[jj]=i  ;this deals with fact that a l-o-s can cross an arm several times (some armas make several turns)
        jj=jj+1
     endif
   endfor
   ;if dlong changes  sign, the arm is crossing the line-of-sight.
   ; we keep in memory ncross, the number of the element of the arm  array where arm is crossed by l-o-s

   ;ncr refers to the points along the arms; kbrac to the points along the line-of-sight
ncr = ncross[0]
if ncr gt 1 then begin
dsun=sqrt(arm[ncr,0]^2 + (arm[ncr,1]-r0)^2)
karm=fix(dsun/step)
 ;karm is the index of the point of the l-o-s where it crosses the arm, will be center of Gaussian distribution

   ;** calculation of angle between arm and l-o-s, in order to obtain the width of the Gaussian
   ;density distribution
      alph=atan((arm[ncr+1,1]-arm[ncr-1,1])/(arm[ncr+1,0]-arm[ncr-1,0]))
   ;delta Y/delta X gives the angle of arm with respect to horizontal axis
      beta=glong-!pi/2
      factor=1./abs(sin(glong+alph))
    if factor ge 10 then factor=10.
    ;this limits the width of the Gaussian to a factor 10 larger than when l-o-s crosses perpendicular

    sig2=2*(factor*sig/.05)^2

    for k=0,nstep-1 do begin
      karmdelt=abs(k-karm)
      if karmdelt gt 15 then begin
        dens[k]=0.
        endif else begin
       dens[k]=denmax*exp( -((k-karm)^2)/sig2) ;fills dens array with Gaussian density distribution
       endelse
    endfor
    endif

;let us consider the possibility that a l-o-s crosses 2 times the same arm, all the same as above
   ncr=ncross[1]
    if ncr gt 1 then begin
    dsun=sqrt(arm[ncr,0]^2 + (arm[ncr,1]-r0)^2)
    karm=fix(dsun/step)

      alph=atan((arm[ncr+1,1]-arm[ncr-1,1])/(arm[ncr+1,0]-arm[ncr-1,0]))
   ;delta Y/delta X gives the angle of arm with respect to horizontal axis
      beta=glong-!pi/2
      gam=abs(glong+alph)
    if gam ge .1 then factor=1./sin(gam) else factor=10.
    ;this limits the width of the Gaussian to a factor 10 larger than when l-o-s crosses perpendicular

    sig2=2*(factor*sig/.05)^2

    for k=0,nstep-1 do begin
      karmdelt=abs(k-karm)
      if karmdelt gt 15 then begin
        dens[k]=0.
        endif else begin
       dens[k]=denmax*exp( -((k-karm)^2)/sig2) ;;k can have any avalue
       endelse
    endfor
    endif

RETURN, dens           ;returns the array with density along l-o-s to main program
END


