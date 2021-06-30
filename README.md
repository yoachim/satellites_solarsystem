# satellites_solarsystem
Checking how mega constellations impact solar system discoveries


Using some code from earlier repos, but whatever

grabbing data files from:
astro-lsst-01.astro.washington.edu:"/local/lsst/opsim/vatiras/"




## Results

If I run on twi_neo_pattern1_v1.7_10yrs veritas, looking at just year 1 with a constellation of 47,000 satellites:

5 arcsec tolerance:
tol_5 0.6136957658128593

15 arcsec tolerance:
tol_15 0.23183481442760062

I think this should be about the same if I do any of the other pattterns, right?
trying twi_neo_pattern3_v1.7_10yrs for year 1:

tol_5 0.6118523115110146
tol_15 0.24790567793980764

Indeed, that's a good match!

