Trying again with pha population

Need to install pycraf on hyak for this


sed -i ''  's|["'\'']||g' survived_merged.txt

to remove all the quotes in a file


run_moving_calc --characterization inner --obsFile survived_merged.txt --opsimDb twi_neo_pattern3_v1.7_10yrs.db --orbitFile /Users/yoachim/rubin_sim_data/orbits/granvik_pha_5k.txt --outDir twi_neo_pattern3_v1.7_10yrs_ss --opsimRun twi_neo_pattern3_v1.7_10yrs --hMin 16.00 --hMax 22.20 --hStep 0.20 --metadata NEO  ; run_moving_fractions --workDir twi_neo_pattern3_v1.7_10yrs_ss --metadata NEO --hMark 22


-----------

Looking at Cumulative completeness at H=22. The full run recovers 55%, the starlink drops to 43% (3 pairs in 12 nights).

original had 2317824 observations, post starlink had 1677992. So we only have 72% of the observations from before--which translates to 78% of the orbits recovered from before.



