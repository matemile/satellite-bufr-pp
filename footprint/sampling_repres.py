#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate  # type: ignore
from scipy.interpolate import griddata
from mpl_toolkits.basemap import Basemap
from eccodes import (
    codes_bufr_new_from_file,
    codes_set,
    codes_get_array,
    codes_release,
)

"""
****************************
Footprint operator functions
****************************
"""

def instrument_geom(
        scan: int,
        sens: int,
):
    """
    Determine instrument geometry necessary for
    the footprint or IFOV size calculation.
    For the time being, AMSU-A, MHS, and IASI sensors
    are set.

    :param scan: Scanposition available from L1b BUFR.
    :param sens: Sensor number (taken from ARPEGE/IFS variational code

    :return: step:   sampling angular interval
             ifov:   IFOV width
             offset: half of the step interval
             nbp:    number of pixels (in half of the full cross scan)
             scr:    modified scan position for the footprint comp.
             denomi: splitter parameter for the sampling resolution
             nprofs: number of footprint operator points
    """

    if sens == 15: # MHS
        step = 1.1111*(np.pi/180.0)
        ifov = 1.1*(np.pi/180.0)
        offset = step/2.0
        nbp = 45
        denomi = 9
        nprofs = 45
        movement_add = (step/denomi)*(np.pi/180.0)
        if scan < 46:
            scr = 46-scan
        else:
            scr = scan -45
    elif sens == 3: #Amsua
        step = 3.3333*(np.pi/180.0)
        ifov = 3.3*(np.pi/180.0)
        offset = step/2.0
        nbp = 15
        denomi = 11
        nprofs = 77
        movement_add = ((step*(180.0/np.pi))/denomi)*(np.pi/180.0)
        if scan < 16:
            scr = 16-scan
        else:
            scr = scan -15
    elif sens == 16: #IASI
        step = 1.65*(np.pi/180.0)
        ifov = 0.84*(np.pi/180.0)
        offset = step/2.0
        nbp = 30
        scr = scan/4.0
    else:
        print('No such sensor definition->exit')
        return

    return step, ifov, offset, nbp, scr, movement_add, denomi, nprofs

def rotation_fov(
        ilat: np.array,
        ilon: np.array,
        aa: float,
        nprofs: int
):
    """
    Rotate geolocations (determined by footprint function)
    by the azimuth angle.

    :param ilat: latitudes of footprint operator points (FOPs).
    :param ilon: longitude of FOPs.
    :param aa: Azimuth angle.
    :param nprofs: Number of FOPs

    :return: Rotated latitudes and longitudes of FOPs.
    """

    #Azimuth
    aa = aa-90.0

    #Rotational matrix
    rot = np.array([[np.cos(aa*np.pi/180.0) , -np.sin(aa*np.pi/180.0)], \
            [np.sin(aa*np.pi/180.0) , np.cos(aa*np.pi/180.0)]])

    #Arrays of rotated geolocations
    rlat = np.zeros(nprofs)
    rlon = np.zeros(nprofs)

    #Distance matrices
    Ell = np.zeros((2,nprofs))
    Ell_rot = np.zeros((2,nprofs))

    #Fill value in the distance matrix
    for kk in range(nprofs):
        Ell[0,kk] = (ilat[kk]*180/np.pi)-(ilat[0]*180/np.pi)
        Ell[1,kk] = (ilon[kk]*180/np.pi)-(ilon[0]*180/np.pi)

    #Do the rotation
    Ell_rot = np.matmul(rot,Ell)

    #Get the rotated geolocations
    for kk in range(nprofs):
        rlat[kk] = Ell_rot[0,kk]+(ilat[0]*180/np.pi)
        rlon[kk] = (Ell_rot[1,kk]/np.cos(ilat[0]))+(ilon[0]*180/np.pi)

    return rlat, rlon

def footprint(
        lat: float,
        lon: float,
        aa2s: float,
        scp: int,
        zsat: float,
        sens: int,
        ra: float = 6378,
        demoni: int = 0,
        nprofs: int = 0,
        movement: float = 0.0,
        gamma: float = 0.0,
        gamma0: float = 0.0,
        gamma1: float = 0.0,
        gamma2: float = 0.0,
        km: float = 0.0,
        kmnorth: float = 0.0,
        path: float = 0.0,
        alpha: float = 0
):
    """
    Footprint operator points

    :param lat: BUFR observation location (latitude).
    :param lon: BUFR observation location (longitude).
    :param aa2s: BUFR Azimuth angle.
    :param scp: BUFR scan position.
    :param zsat: BUFR satellite altitude.
    :param sens: Sensor id (ARPEGE/IFS code)
    :param ra: Earth's radius in km
    :param denomi: splitter parameter for the sampling resolution
    :param nprofs: number of footprint operator points
    :param alpha: viewing or scanning angles
    :param gamma*: scanning angles on the Earth's frame
    :param km: distance (in km) of IFOV major axis
    :param kmnorth: distance (in km) if IFOV minor axis
    :param path: slant-path distance

    :return: Rotated latitudes and longitudes of FOPs.
    """

    # Collect details of instrument geometry
    step, ifov, offset, nbp, scanpos, movement, denomi, nprofs = instrument_geom(scp, sens)

    #Compute the angles of the given IFOV ellipse
    #Viewing angles of IFOV ellipse cross-track vertices (1-start, 2-end)
    alpha = offset + (scanpos-1)*step - ifov/2.
    gamma1 = np.arcsin((ra+zsat)/ra * np.sin(alpha)) - alpha
    alpha = offset + (scanpos-1)*step + ifov/2.
    gamma2 = np.arcsin((ra+zsat)/ra * np.sin(alpha)) - alpha
    #Bore-sight viewing angles
    alpha = offset + (scanpos-1)*step
    gamma = np.arcsin((ra+zsat)/ra * np.sin(alpha)) - alpha
    # Distance of IFOV major axis on the ground
    km = (gamma2 - gamma1)*ra
    path = ra*np.sin(gamma)/np.sin(alpha)
    # Distance of IFOV minor axis on the ground
    kmnorth = path*ifov
    km *= 1.25 #1.0+(step/11.0) Experimental threshold...

    print(f"Footprint-op-Along-track: {kmnorth:.2f} {'km'}")
    print(f"Footprint-op-Across-track: {km:.2f} {'km'}")
    print(f"Footprint-op-Area: {(km/2.0)*(kmnorth/2.0)*np.pi:.2f} {'km2'}")

    plat = np.zeros(nprofs)
    plon = np.zeros(nprofs)
    
    #Bore-sight position
    plat[0] = lat*np.pi/180.0
    plon[0] = lon*np.pi/180.0
    
    ii = 0

    #Angles to be used for sampling
    circd = [0.0, 90.0, 45, 22.5, 67.5, 11.25, 33.75, 56.25, 78.75]

    for i in range(1,int(denomi/2.0)+1):

        #Sampling resolution steps
        aa = i*kmnorth/denomi
        bb = i*km/denomi #+ (movement*ra)
        #print('aa',aa)
        #print('i-bb',i,bb,movement,(movement*ra))

        #Thinning the sampling near the bore-sight
        if i == 1: limd = 2
        if i == 2: limd = 3
        if i == 3 or i == 4 or i == 5 or i == 6: limd = 5
        if i >= 7: limd = 9

        for iid in range(0,limd):

            ddeg = circd[iid]
            ddeg *= np.pi/180

            plat[ii+1] = plat[0] + aa*(np.sin(ddeg))/ra
            plon[ii+1] = plon[0] + bb*(np.cos(ddeg))/ra

            plat[ii+2] = plat[0] - aa*(np.sin(ddeg))/ra
            plon[ii+2] = plon[0] - bb*(np.cos(ddeg))/ra

            if circd[iid] == 0.0 or circd[iid] == 90.0:
                ii += 2
                continue

            plat[ii+3] = plat[0] - aa*(np.sin(ddeg))/ra
            plon[ii+3] = plon[0] + bb*(np.cos(ddeg))/ra

            plat[ii+4] = plat[0] + aa*(np.sin(ddeg))/ra
            plon[ii+4] = plon[0] - bb*(np.cos(ddeg))/ra

            ii += 4

    return rotation_fov(plat, plon, aa2s, nprofs)

def read_single_obs_metadata(bufr_file):
    """
    Read metadata from a single‑observation BUFR file.

    Returns:
        lat, lon, zenith, azimuth, scanp, zsat  (all floats/ints)
    """
    lat_key     = "latitude"
    lon_key     = "longitude"
    zenitha_key = "satelliteZenithAngle"
    azimutha_key = "bearingOrAzimuth"
    scanp_key   = "fieldOfViewNumber"
    zsat_key    = "heightOfStation"

    with open(bufr_file, "rb") as f:
        h = codes_bufr_new_from_file(f)
        if h is None:
            raise RuntimeError("No BUFR message found")

        codes_set(h, "unpack", 1)

        # All are 1‑element arrays in your single‑subset file → take [0]
        lat     = float(codes_get_array(h, lat_key)[0])
        lon     = float(codes_get_array(h, lon_key)[0])
        zenith  = float(codes_get_array(h, zenitha_key)[0])
        azimuth = float(codes_get_array(h, azimutha_key)[0])
        scanp   = float(codes_get_array(h, scanp_key)[0])
        zsat    = float(codes_get_array(h, zsat_key)[0])

        codes_release(h)

    return lat, lon, zenith, azimuth, scanp, zsat

def plot_on_map(lat, lon, flat, flon, sensname):
    """
    Plot antenna pattern and footprint representation on a Basemap 
    with an automatically chosen domain around the observation boresight. 
    """

    #fig, ax = plt.subplots()
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    
    #Use a different projection
    map = Basemap(projection='stere', #Single IFOV scl9
            lon_0 = lon,
            lat_0 = lat,
            width = 200000,
            height = 200000,
            resolution='i')
    
    # draw parallels
    parallels = np.arange(0.,90,1.0)
    map.drawparallels(parallels,labels=[1,1,0,0],fontsize=8)#, rotation=90)
    # draw meridians
    meridians = np.arange(0,360.,2.0)
    map.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)#, rotation=90)
    # Coastlines
    map.drawcoastlines(linewidth=1.0)

    cmap = plt.cm.get_cmap("Greys",10)
    #map.scatter(w[:, 1], w[:, 0], c=w[:, 2], latlon=True, s=50, cmap=cmap)
    map.scatter(flon, flat, latlon=True, s=30, color='r', marker='.')
    plt.title("Footprint and its operator points on the ground \n Sensor: "+sensname+", Lat:"+str(lat)+", Lon:"+str(lon)+"")
    plt.savefig(""+sensname+"_antenna_sampling.png")

#
# Main
#

if len(sys.argv) < 4:
    print(f"Usage: {sys.argv[0]} BUFR_FILE SENSOR SENSORNAME", file=sys.stderr)
    print(
        f"{sys.argv[0]} BUFR_FILE SENSOR SENSORNAME",
        file=sys.stderr,
    )
    sys.exit(1)

bufr_file = sys.argv[1]
sens = int(sys.argv[2])
sensname = str(sys.argv[3])

lat, lon, zenith, azimuth, scanp, zsat = read_single_obs_metadata("../singleobs/single_subset_"+sensname+".bufr")

print('lat, lon, zenith, azimuth, scanp, zsat',lat, lon, zenith, azimuth, scanp, zsat)

#Compute footprint operator points on the ground
flat, flon = footprint(lat, lon, azimuth, scanp, zsat/1e3, sens)

print('flat, flon',flat, flon)

# Plot antenna pattern and footprint operator points on map
plot_on_map(lat, lon, flat, flon, sensname)
