## NMSYN

#### To Do list: 
- [ ] Test gravity strain traces in surface directions
- [ ] Add a flag for rotation ZNE or RTZ output instead of manual

#### Updates to the code: 
- Most recent version allows for gravity strain tensor computation. I note here that for the tensor to be unitless, it is typically defined as $$\mathbf{H}(\mathbf{r}, t) = \int^t_0 \int^{t'}_0 \nabla \nabla \phi(\mathbf{r},t'') \: \mathrm{d}t''\mathrm{d}t. $$ Here, the code outputs only the components of $$ \nabla \nabla \phi(\mathbf{r},t'')\quad, $$ that is you will need to do the time integration separately. 
- The code now contains a flag ```grav_inside``` defined in ```spherical.h```. This allows you to switch the definition of the gravitational acceleration and strain dependent on whether you are want the value just above/below the surface (since these values are discontinuous across the surface). 

### Introduction

Originally Developed by: Jeroen Tromp

Modified by: \
03/21/02:  Qinya Liu (liuqy@physics.utoronto.ca) \
02/11/22:  Hom Nath Gharti (hng1@queensu.ca) \
05/20/24:  Will Eaton (weaton@princeton.edu)\ 

```xnmsyn``` calculates a normal mode synthetics using normal mode summation for given ```CMTSOLUTION``` and ```STATIONS``` files. 
It relies on binary databases of modes computed using [MINEOS](https://github.com/geodynamics/mineos)
For editors of the code, it is recommended to use Dahlen and Tromp (1998) as a reference for equations used. 

#### 1. INPUT PARAMTERS:

EXAMPLE--- 
```
./database/prem_noocean_model
./database/toroidal_prem_noocean.bin
./database/radial_prem_noocean.bin
./database/spheroidal_prem_noocean.bin
./input/CMTSOLUTION
./input/STATIONS
17.0
1 12000
```

are respectively: model deck file, toroidal catalogue, radial catalogue,
and spheroidal catalogue, CMTSOLUTION file name, STATIONS file name, shortest period 
in summation, dt in output time series, and total length of record in 
seconds.

example execution: ```xnmsyn input/prem_noocean.in```

#### 2. SETUP

  Make sure you have output/ folder in the current directory.

file STATIONS has the format of
```
2
KMI  IC   25.1233   102.7400 1975.0    0.0
PAS  TS   34.1490  -118.1720    0.0    0.0
```
first line -- number of stations
second - last line -- station_name, network_code stalat, stalon, staele, burial


CMTSOLUTION files are the Harvard format:
```
PDE  1994  6  9  0 33 45.40 -13.8200  -67.2500 647.1 6.9 6.8      NORTHERN BOLIVIA
event name:    M060994A
time shift:      0.0000
half duration:  30.0000
latitude:      -13.8200
longitude:     -67.2500
depth:         647.1000
Mrr:      -0.759000E+28
Mtt:       0.775000E+28
Mpp:      -0.160000E+27
Mrt:      -0.250300E+29
Mrp:       0.420000E+27
Mtp:      -0.248000E+28
```
#### 3. OUTPUT 

OUTPUT will be in output/ with names like:
```PAS.TS.LHE.nmsyn  PAS.TS.LHN.nmsyn  PAS.TS.LHZ.nmsyn```
and the old version of synthetics with the same name will
be overwritten, so be sure to rename the old ones if they
are still useful.

You can type 
```xmgr -graph 0 PAS.TS.LHE.nmsyn -graph 1 PAS.TS.LHN.nmsyn -graph 2 PAS.TS.LHZ.nmsyn -rows 3``` to see preliminary result.


### REMARKS:

1. If you want to make synthetics with another mode catalogue, you might
want to change the parameters in minos.h, including NR and NL. Make sure
you find out the written format of the catalogue and change correspondingly
in the nmsyn.f.

2. If you want output to be RTZ component instead of NEZ component, you can
change parts of write_record.f, or do the rotation with your own program.

3. A minor problem exists for the output synthetics. The synthetics is not
perfectly zero before P arrivals, and I cannot figure out what causes this.
You are more than welcomed to find out why. ( Jeroen has a program to remove
the offset)

4. There are several mode catalogues available for summation:
```
radial.aniprem808_40s.bin
toroidal.aniprem808_40s.bin
spheroidal.aniprem808_40s.bin
``` 
is a 40s mode catalogue for aniprem808 with top 165 layers stored.

```
radial.aniprem808_8s.bin
toroidal.aniprem808_8s.bin
spheroidal.aniprem808_8s.bin
```
is a 8s mode catalogue for aniprem808 with top 165 layers stored.

Use the 40s catalogue directly instead of using tmin=40 for the 8s
catalogue will save you a lot of time.

5. One can compare the output with the output of QmXD in Sun
workstation, using the following script

```
#!/usr/bin/tcsh -f
setenv MYDBS /net/freawaru/export/home2/jtromp/dbs
/net/freawaru/export/home2/jtromp/bin/QmXD -M \
/net/freawaru/export/home1/jtromp/dbs/catalogues/seigsml_prem_an.750.8s \
-O 8. -P 6000. -s
```
The output generally agrees farely well except the problem I described
earlier.

6. The dt can be chosen as long as 2*dt < Tmin, larger dt will result
in faster computation.

7. The program was originally written by Jeroen, and I only changed the
input and output format.


