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

-------------

ok, now to try running the original and the surviving through MAF. 
Kind of a mess of things since this was run in parallel for some reason.

grabbing the initial orbit file from hyak:/gscratch/astro/lynnej/orbits/vatiras_granvik_10k


The initial commands look something like this. 
run_moving_calc.py --orbitFile /gscratch/astro/lynnej/orbits/vatiras_granvik_10k/vatiras_granvik_10k.txt --opsimRun twi_neo_pattern3_v1.7_10yrs --opsimDb /gscratch/astro/lynnej/fbs_1.7/twi_neo/twi_neo_pattern3_v1.7_10yrs.db --characterization inner --hMin 16.0 --hMax 28.0 --hStep 0.2 --metadata Vatira --nYearsMax 10 --startTime 59853 --outDir vatiras_granvik_10k_2  --obsFile twi_neo_pattern3_v1.7_10yrs__vatiras_granvik_10k_2_obs.txt

This looks like it works:

run_moving_calc --orbitFile orbit_files/vatiras_granvik_10k.txt --opsimRun twi_neo_pattern3_v1.7_10yrs --opsimDb twi_neo_pattern3_v1.7_10yrs.db --characterization inner --hMin 16.0 --hMax 28.0 --hStep 0.2 --metadata Vatira --nYearsMax 10 --startTime 59853 --outDir vatiras_granvik_10k_init  --obsFile init_twi_neo_pattern3_v1.7_10yrs__vatiras_granvik_10k.txt

