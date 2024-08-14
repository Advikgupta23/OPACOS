;glong = 67.00
;glat = -13.00

;glong = 270.0
;glat = 4.0

glong = 270.0
glat = -73.0

; Define the range of distances (in kpc) you want to loop through
start_dist = 0.01  ; Starting distance
end_dist = 9.92   ; Ending distance
step_dist = 0.1   ; Step size

; Calculate the number of steps
nsteps = fix((end_dist - start_dist) / step_dist) + 1

; Create arrays to store the results
distances = make_array(nsteps, /float)
ebvspirals = make_array(nsteps, /float)
avspiral = make_array(nsteps, /float)

; Loop through the distances and call extinspiral
for i = 0, nsteps-1 do begin
    dist = start_dist + i * step_dist
    distances[i] = dist
    extin, glong, glat, dist, ebvspiral, avspiral[i]
    ebvspirals[i] = ebvspiral
    ;print, ebvspiral

endfor

; Print the results
;print, 'EBV Spiral'
for i = 0, nsteps-1 do begin
    print, ebvspirals[i],','
endfor

end

