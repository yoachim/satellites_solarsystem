from utils import Constellation, starlink_constellation
import numpy as np
import glob

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
    result = np.genfromtxt(filename, names=True, skip_header=14)
    return result


def run_check(path='baseline_nexp2_v1.7_10yrs', streak_tol=20., start_dist=10.):
    """
    Paramters
    ---------
    streak_tol : float (20)
        The tolerance to alow objects to be from a streak (arcsec)
    start_dist : float (10)
        How far away to consider a satellite a threat in computing streak (degrees)
    """

    streak_tol = np.radians(streak_tol/3600.)  # to radians
    start_dist = np.radians(10)

    files = glob.glob(path+'/*.txt')
    arrays = [read_ss_objs(filename) for filename in files]
    ss_data = np.hstack(arrays)
    # XXX---starting with small constellation for speed to start
    sats = starlink_constellation(fivek=True, alt_limit=20.)#(supersize=True)
    const = Constellation(sats)


    # I think I want to output mjd, ra, dec for the satellites that get hit.
    for mjd_start in np.unique(ss_data['time']):
        const.update_mjd(mjd_start)
        alt_rad_start = const.altitudes_rad[const.above_alt_limit] + 0
        az_rad_start = const.azimuth_rad[const.above_alt_limit] + 0
        obj_indxs = np.where(ss_data['time'] == mjd_start)
        for obj_indx in obj_indxs:
            pass

if __name__ == '__main__':
    run_check()


    


