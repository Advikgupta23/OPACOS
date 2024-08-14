Function spiral2, r0,fase,pitch0,dpitch,npont
; this function creates a 2-dimensional array containing NPONT lines, one for each point of the arm
; the columns are X,Y,R,THETA coordinates of each point (all referred to the Galactic center)
; parameters are: r0= initial radius of the arm , fase= initial azimuthal angle (direction of Sun is fase=0)
;pitch0= initial pitch angle of the arms (angle with respect to circles). for most arms this is constant,
; but we introduced the possibility of variable (steadily increasing or decreasing) pitch angle.
;dpitch is the rate of variation

pitch = make_array(npont)
r=make_array(npont)
teta=make_array(npont)
x=make_array(npont)
y=make_array(npont)

conv=!pi/180
cte=tan(pitch0*conv)
tet=findgen(npont)*conv ;creates an array with theta for each point of the arm, they are spaced one degree
faserad=fase*conv
teta=tet+faserad
r[0]=r0
pitch[0]=pitch0
for j=1,npont-1 do begin
     pitch[j]=pitch0+j*dpitch      ;the pitch angle can vary!
     r[j]=r[j-1] *(1.+ tan(pitch[j-1]*conv)*conv)
   endfor
;above, for each radius, the next radius is the old one plus the increase in radius,
;which is equal to the delta angle (one degree=conv)times the tangent of the pitch angle
;r=r0*EXP(teta*cte) is the usual equation of an arm; we do not use it to allow variable pitch
x=-r*sin(teta)
y=r*cos(teta)
spi=[[x],[y], [r], [teta]]
return,spi   ;returns the arm array to the main program
end


