#!/usr/bin/env python

from utils import Constellation, starlink_constellation
import numpy as np
import glob
from rubin_sim.utils import _approx_RaDec2AltAz, Site
import sys
import pandas as pd
import argparse

# Bearing and track length from: 
# https://stackoverflow.com/questions/32771458/distance-from-lat-lng-point-to-minor-arc-segment


def _angularSeparation(long1, lat1, long2, lat2):
    """
    """
    t1 = np.sin(lat2/2.0 - lat1/2.0)**2
    t2 = np.cos(lat1)*np.cos(lat2)*np.sin(long2/2.0 - long1/2.0)**2
    _sum = t1 + t2

    if np.size(_sum) == 1:
        if _sum < 0.0:
            _sum = 0.0
    else:
        _sum = np.where(_sum < 0.0, 0.0, _sum)

    return 2.0*np.arcsin(np.sqrt(_sum))


def bear(latA, lonA, latB, lonB):
    """All radians
    """
    # BEAR Finds the bearing from one lat / lon point to another.
    result = np.arctan2(np.sin(lonB - lonA) * np.cos(latB),
                        np.cos(latA) * np.sin(latB) - np.sin(latA) * np.cos(latB) * np.cos(lonB - lonA)
                        )

    return result


def pointToLineDistance(lon1, lat1, lon2, lat2, lon3, lat3):
    """All radians
    points 1 and 2 define an arc segment,
    this finds the distance of point 3 to the arc segment. 
    """

    result = lon1*0
    needed = np.ones(result.size, dtype=bool)

    bear12 = bear(lat1, lon1, lat2, lon2)
    bear13 = bear(lat1, lon1, lat3, lon3)
    dis13 = _angularSeparation(lon1, lat1, lon3, lat3)

    # Is relative bearing obtuse?
    diff = np.abs(bear13 - bear12)
    if np.size(diff) == 1:
        if diff > np.pi:
            diff = 2*np.pi - diff
        if diff > (np.pi / 2):
            return dis13
    else:
        solved = np.where(diff > (np.pi / 2))[0]
        result[solved] = dis13[solved]
        needed[solved] = 0
    
    # Find the cross-track distance.
    dxt = np.arcsin(np.sin(dis13) * np.sin(bear13 - bear12))

    # Is p4 beyond the arc?
    dis12 = _angularSeparation(lon1, lat1, lon2, lat2)
    dis14 = np.arccos(np.cos(dis13) / np.cos(dxt))
    if np.size(dis14) == 1:
        if dis14 > dis12:
            return _angularSeparation(lon2, lat2, lon3, lat3)
    else:
        solved = np.where(dis14 > dis12)[0]
        result[solved] = _angularSeparation(lon2[solved], lat2[solved], lon3[solved], lat3[solved])

    if np.size(lon1) == 1:
        return np.abs(dxt)
    else:
        result[needed] = np.abs(dxt[needed])
        return result


def read_ss_objs(filename):
    dt = [('objId', '<U8'), ('time', '<f8'), ('ra', '<f8'), ('dec', '<f8'),
    ('dradt', '<f8'), ('ddecdt', '<f8'), ('phase', '<f8'), ('solarelon', '<f8'),
    ('helio_dist', '<f8'), ('geo_dist', '<f8'), ('magV', '<f8'), ('trueAnomaly', '<f8'),
    ('velocity', '<f8'), ('fieldDec', '<f8'), ('fieldRA', '<f8'), ('filter', '<U1'),
    ('fiveSigmaDepth', '<f8'), ('night', '<f8'), ('observationStartMJD', '<f8'),
    ('rotSkyPos', '<f8'), ('seeingFwhmEff', '<f8'), ('seeingFwhmGeom', '<f8'), ('solarElong', '<f8'),
    ('visitExposureTime', '<f8'), ('dmagColor', '<f8'), ('dmagTrail', '<f8'), ('dmagDetect', '<f8')]
    result = np.genfromtxt(filename, skip_header=0, dtype=dt)
    return result


def run_check(infile, streak_tol1=5., streak_tol2=15., year1=True):
    """
    Paramters
    ---------
    streak_tol : float (20)
        The tolerance to alow objects to be from a streak (arcsec)
    start_dist : float (10)
        How far away to consider a satellite a threat in computing streak (degrees)
    """
    outfile = 'survived_' + infile

    site = Site(name='LSST')
    streak_tol1 = np.radians(streak_tol1/3600.)  # to radians
    streak_tol2 = np.radians(streak_tol2/3600.)
    print('reading initial')
    ss_data = read_ss_objs(infile) 

    ss_data.sort(order='time')
    
    # Let's crop down to first year:
    if year1:
        good = np.where(ss_data['time'] <= np.min((ss_data['time']+365.25)))[0]
        ss_data = ss_data[good]

    # XXX---starting with small constellation for speed to start
    sats = starlink_constellation(supersize=True)
    const = Constellation(sats, alt_limit=20.)
    const.advance_epoch(advance=100)

    # array to hold if something is too close to a detection
    ss_obj_visible1 = np.ones(ss_data.size, dtype=bool)
    ss_obj_visible2 = np.ones(ss_data.size, dtype=bool)

    _ones = np.ones(ss_data.size)

    ss_alt, ss_az = _approx_RaDec2AltAz(np.radians(ss_data['fieldRA']), np.radians(ss_data['fieldDec']),
                                        site.latitude_rad*_ones, site.longitude_rad*_ones, ss_data['time'])

    # I think I want to output mjd, ra, dec for the satellites that get hit.
    #n_max = np.unique(ss_data['time']).size
    for i, mjd_start in enumerate(np.unique(ss_data['time'])):
        const.update_mjd(mjd_start)
        alt_rad_start = const.altitudes_rad + 0
        az_rad_start = const.azimuth_rad + 0
        obj_indxs = np.where(ss_data['time'] == mjd_start)[0]
        # Assume they all have the same visit Exposure time
        const.update_mjd(mjd_start + ss_data[obj_indxs[0]]['visitExposureTime'])
        alt_rad_end = const.altitudes_rad
        az_rad_end = const.azimuth_rad
        for obj_indx in obj_indxs:
            distances = pointToLineDistance(az_rad_start, alt_rad_start, az_rad_end, alt_rad_end,
                                            ss_az[obj_indx]*np.ones(az_rad_start.size),
                                            ss_alt[obj_indx]*np.ones(az_rad_start.size))
            if np.min(distances) <= streak_tol1:
                ss_obj_visible1[obj_indx] = False
            if np.min(distances) <= streak_tol2:
                ss_obj_visible2[obj_indx] = False
        #progress = i/n_max*100
        #text = "\rprogress = %.1f%% (%i of %i)" % (progress, i, n_max)
        #sys.stdout.write(text)
        #sys.stdout.flush()
    result = {'tol_%i' % (np.degrees(streak_tol1)*3600): ss_obj_visible1,
              'tol_%i' % (np.degrees(streak_tol2)*3600): ss_obj_visible2}

    ss_data = ss_data[ss_obj_visible1]
    _tempdf = pd.DataFrame(ss_data)
    _tempdf.to_csv(outfile, index=False, sep=' ')

    return result


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="check which asteroid observations hit a streak")
    parser.add_argument("infile", type=str)
    args = parser.parse_args()

    infile = args.infile
    vis = run_check(infile, year1=False)
    for key in vis:
        print(key, np.mean(vis[key]))

