# SN II-P fitting code
### A global model of the light curves and expansion velocities of Type II-Plateau Supernovae

Code by Ondřej Pejcha (Charles University, Czech Republic) -- e-mail: pejcha@utf.mff.cuni.cz

Minor modifications by Tomás Müller (U. of Southampton, UK) -- e-mail: t.e.muller-bravo@soton.ac.uk

README by Tomás Müller


## Important Note

This code was developed by Ondřej Pejcha (Charles University) and José Luis Prieto (Universidad Diego Portales) code to fit SNe II light curves
and velocities. I made this repository as his Ondřej's website from Princeton University does not longer exists.

I modified Ondřej's version to give output files that are easier to handle. I also added the option to manually 
give an explosion epoch initial estimation and errors constrains. It is easy to do the same for other parameters 
by following the changes I made (check `fit_single.cpp` file). The code now rescales the errors (with <img src="https://latex.codecogs.com/gif.latex?\chi^2"/>) 
before estimating the physical parameters with the analytical relations from Litvinova & Nadëzhin (1985) and Popov (1993). 
These output files have slightly different names from previous versions of the code.

There are comments added with the changes made, but let me know if there is anything not clear.
An `index.html` file is included with information from the Ondřej's website from Princeton.

I created a python script (`show_output.py`) to plot the output fits from the modified code 
(run this from the top directory, `src` by default). Thishelps to have a quick look and see 
if the fits are working correctly. You have to run it from the code's top directory (i.e., 
where the `fit_single.cpp` file is) and it will ask you for the SN name. You can run this by typing:

```
$ python show_output.py
```

## Installation 

To build the package you need to follow these steps:

```
$ cd /<path_to>/cmpfit-1.2
$ make clean
$ make
$ cd ..
$ make fit_single
```

and that is all. You can try it with an example file it comes along the installation package:

```
$ ./fit_single.exe SN2013am.dat
```

The code will ask you to give the SN name. If you input `SN2013am`, a directory with that exact name 
will be created and all the output files of the SN will end there. If it doesn't complain or give you any error, you are ready!
**Note** that in order for the code to work, apart from the photometry, at least one velocity measurement should be included.

## Issues and Contributing

If you would like to contribute or raise issues, you can either 1) submit an [issue on Github](https://github.com/temuller/snii_fitting_code/issues)
(and submit a pull request if you wish) or 2) drop me, Ondřej or both an email.


## Input File 

These are the columns the input files must have: ref dset flt mjd val err. See below a more detailed description of these:

| column  |  dtype   |                    Description                           |
|:-------:|:--------:|:--------------------------------------------------------:|
|  ref 	  |  integer |this does not really do anything, at least that I know    |
|  dset	  |  integer |1 for photometry or 0 for velocities (FeII 5169 only)	|
|  flt 	  |  integer |filter number (see below). If dset==0, set it to 0        |
|  mjd 	  |   float  |modified julian date [MJD]                                |
|  val 	  |   float  |magnitudes (if dset==1) or velocities in m/s (if dset==0) |
|  err 	  |   float  |error of 'val'                                            |


Filter numbers for **flt** column:

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


## License & Attribution

Copyright 2014-2020 Ondřej Pejcha and José Luis Prieto.

The source code is made available under the terms of the MIT license.

If you make use of this code, please cite the paper which is currently on [ADS](https://ui.adsabs.harvard.edu/abs/2015ApJ...799..215P/abstract):

```
@ARTICLE{2015ApJ...799..215P,
       author = {{Pejcha}, Ond{\v{r}}ej and {Prieto}, Jose L.},
        title = "{A Global Model of The Light Curves and Expansion Velocities of Type II-plateau Supernovae}",
      journal = {\apj},
     keywords = {methods: statistical, stars: distances, supernovae: general, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
         year = 2015,
        month = feb,
       volume = {799},
       number = {2},
          eid = {215},
        pages = {215},
          doi = {10.1088/0004-637X/799/2/215},
archivePrefix = {arXiv},
       eprint = {1409.2500},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015ApJ...799..215P},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```


