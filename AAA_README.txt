Code by Ondřej Pejcha (Charles University, Czech Republic)

# March 2020
Minor modifications by Tomás Müller (U. of Southampton, UK)
e-mail: t.e.muller-bravo@soton.ac.uk

README by Tomás Müller

####################
## Important Note ##
####################
This is Ondrej Pejcha's code to fit SNe II light curves. I made this as his Princeton University's website 
does not longer exists.

I modified Ondrej's version to give output files that are easier to handle. I also added the option to manually 
give an explosion epoch initial estimation and errors constrains. It is easy to do the same for other parameters 
by following the changes I made (check "fit_single.cpp" file). The code now rescales the errors (with Chi^2) 
before estimating the physical parameters with LN85 and POPOV93. These output files have slightly different 
names from previous versions.

There are comments added with the changes made, but let me know if there is anything not clear.
An "index.html" file is included with information from the Ondrej's website from Princeton.

I created a python script (show_output.py) to plot the output fits from the modified code. This
helps to have a quick look and see if the fits are working correctly. You have to run it from the
code's top directory (where the fit_single.cpp file is).

####################
### Installation ###
####################
To build the package you need follow these steps:

$cd /<path_to>/cmpfit-1.2
$make clean
$make
$cd ..
$make fit_single

and that is all. You can try it with an example file it comes along the installation package:

$./fit_single.exe SN2013am.dat (or sn2013am.proc from the previous version)

If it doesn't complain or give you any error, you are ready!

####################
#### Input File ####
####################

These are the columns the input files must have: ref dset flt mjd val err. See below a more detailed description of these:

ref: integer - this does not really do anything, at least that I know.
dset: integer (1 or 0) -  1 for photometry and 0 for velocities (FeII 5169 only).
flt: integer - filter number (see below). If dset==0, set it to 0.
mjd: float - modified julian date [MJD] (I believe it also works with julian dates [JD]).
val: float - magnitudes (if dset==1) or velocities in m/s (if dset==0).
err: float - error of 'val'.


############################################
Filter numbers for flt column:

U ----> 0

B ----> 1

V ----> 2

R ----> 3

I ----> 4

J ----> 5

H ----> 6

K ----> 7

g ----> 9 #SDSS DR7 SDSS = AB - 0.02

r ----> 10 #SDSS DR7 SDSS = AB - 0.02

i ----> 11 #SDSS DR7 SDSS = AB - 0.02

z ----> 12 #SDSS DR7 SDSS = AB - 0.02

uvw2 ----> 13 #Swift uvw2 @ 1928A

uvm2 ----> 14 #Swift uvm2 @ 2246A

uvw1 ----> 15 #Swift uvw1 @ 2600A

u ----> 16 #Swift u    @ 3465A

b ----> 17 #Swift b    @ 4392A

v ----> 18 #Swift v    @ 5468A

Z ----> 19 # Hamuy et al. (2001) - on 1999em

Y ----> 20 # in the YJHK series?
############################################

