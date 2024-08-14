EXTIN.PRO corresponds to the axis-symmetric model of Amores & Lepine
(AJ, 130, 679).
It is written in IDL language and does not require external functions.
Is computes the extintion in the visible (Av) along a path from the
Sun to any point in the Galaxy, specified by galactic coordinates and
distance. It asks for galactic longitude and latitude,(in degrees)
and distance (in kpc).

EXTINSPIRAL.PRO corresponds to the Spiral (S) model of the same paper.
It asks for the same input parameters of EXTIN.PRO and gives the same type
of result. It requires the functions SPIRAL2.PRO and ARMDENS.PRO 
to work. These files must be in the same directory and compiled together.

Additional explanations on the steps performed by the programs have
been introduced as comments in the program themselves.

