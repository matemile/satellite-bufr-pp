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
Basic, stand-alone, functions related to 2D antenna responses.

Basic means that the antenna pattern is narrow and circular at the
antenna, and that a spherical planet can be assumed. This allows to
neglect second order effects.

These coordinates are used
-----
za = zenith angle
aa = azimuth angle, counted clockwise, with zero at North or flight direction
ac = across-track distance
al = along-track distance
lat = latitude
lon = longitude

For za and aa degrees are used (i.e. not radians).

Please note that za, aa, ac and al can be both absolute positions and
deviations to the bore-sight.

Cross-track scanning
-----
The scanning is assumed to be perpendicular to the flight direction.
Antenna patterns expressed in (za,aa) or (ac,al) are valid for the
right side. That is is at absolute aa of 90, going towards positive
za. This means that za and ac coincide, while aa and al are parallel
but go in opposite directions. The movement during integration goves a
broadening in za.

Conical scanning
----
Antenna patterns expressed in (za,aa) or (ac,al) match looking
straight forward (absolute aa of 0). Here za coincides with al, and aa
coincides with ac. The movement during integration goves a broadening
in aa.
"""

def fwhm_by_interp(
    x: np.ndarray,
    y: np.ndarray,
    R: np.ndarray
    ): #-> tuple[float, float]:
    """
    Determines FWHM by interpolation.

    :param x: x-grid, either za or ac.
    :param y: y-grid, either aa or al
    :param R: Input response.

    :return: Tuple with FWHM in x- and y-direction
    """

    # Find index of bore-sight
    x0 = np.where(x == 0)[0]
    if len(x0) != 1:
        raise ValueError("Array x must contain 0.")
    y0 = np.where(y == 0)[0]
    if len(y0) != 1:
        raise ValueError("Array y must contain 0.")
    x0 = int(x0[0])
    y0 = int(y0[0])

    hv = 0.5 * R[x0, y0]

    # x
    i2 = x0 + 1
    F = interpolate.interp1d(R[0:i2, y0], x[0:i2])
    fwhm_x = -F(hv)
    i2 = len(x) - 1
    F = interpolate.interp1d(R[x0:i2, y0], x[x0:i2])
    fwhm_x += F(hv)

    # y
    i2 = y0 + 1
    F = interpolate.interp1d(R[x0, 0:i2], y[0:i2])
    fwhm_y = -F(hv)
    i2 = len(y) - 1
    F = interpolate.interp1d(R[x0, y0:i2], y[y0:i2])
    fwhm_y += F(hv)

    return fwhm_x, fwhm_y


def fwhm2si(
    fwhm: float
) -> float:
    """
    Calculates the Gaussian standard deviation matching a FWHM.

    :param fwhm: The full width at hald max (FWHM).

    :return: The matching 1 std dev.
    """

    return fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def gaussian2d(
    fwhm_x: float,
    fwhm_y: float,
    width_x: float = -1.0,
    width_y: float = -1.0,
    dx: float = -1.0,
    dy: float = -1.0
    ): #-> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generates a 2D Gaussian response.

    The output is a bivariate normal distribution, covering a
    rectangular area.

    Generated grids have spacing specified, while the width can be
    somewhat larger than selected value.

    Units of grids and response returned follow from unit of input.

    :param fwhm_x:  The FWHM in x-direction, matching either za or ac.
    :param fwhm_x:  The FWHM in y-direction, matching either za or ac.
    :param width_x: The minimum width of x-grid. If <= 0, set to be
                    three times the fwhm_x. The later is default.
    :param width_y: The minimum width of y-grid. If <= 0, set to be
                    three times the  fwhm_y. The later is default.
    :param dx:      Spacing of x-grid. If <= 0, set to be 0.05 of
                    fwhm_x. The later is default.
    :param dy:      Spacing of y-grid. If <= 0, set to be 0.05 of
                    fwhm_y. The later is default.

    :return: Tuple with x-grid, y_grid, and response.
    """

    # Use default values?
    if dx <= 0:
        dx = 0.05 * fwhm_x
    if dy <= 0:
        dy = 0.05 * fwhm_y
    if width_x <= 0:
        width_x = 3.0 * fwhm_x
    if width_y <= 0:
        width_y = 3.0 * fwhm_y

    # Create grids
    n = int(np.ceil(width_x / (2.0 * dx)))
    hwidth = n * dx
    grid_x = np.linspace(-hwidth, hwidth, 2 * n + 1)
    n = int(np.rint(np.ceil(width_y / (2.0 * dy))))
    hwidth = n * dy
    grid_y = np.linspace(-hwidth, hwidth, 2*n + 1)

    # Calculate std devs matching FWHMs
    si_x = fwhm2si(fwhm_x)
    si_y = fwhm2si(fwhm_y)

    # Init response array with normalisation term
    nx = len(grid_x)
    ny = len(grid_y)
    R = np.full([nx, ny], 1 / (2.0 * np.pi * si_x * si_y))

    # Add exponential term
    for x in range(nx):
        tmp = np.power(grid_x[x] / si_x, 2)
        for y in range(ny):
            R[x, y] *= np.exp(-0.5 * (tmp + np.power(grid_y[y] / si_y, 2)))

    return grid_x, grid_y, R


def geopos_response(
    ac: np.ndarray,
    al: np.ndarray,
    Rg: np.ndarray,
    lat0: float,
    lon0: float,
    aa2sat: float,
    is_conical: bool = False,
    re: float = 6378.1e3,
    scn: int = 30
) -> np.ndarray:
    """
    Geolocates a response around given bore-sight and angle.

    The input response is assumed to be m,apped to distances at
    ground level.

    :param ac:         Across-track grid.
    :param al:         Along-track grid.
    :param R:          Response mapped to ground distances
    :param lat0:       Latitude of bore-sight
    :param lon0:       Longitude of bore-sight
    :param aa2sat:     Azimuth angle, from ground to satellite [-180, 360].
    :param is_conical: Cross-track (False, defualt) or concical (True) scanner.
    :param re:         Planet radius.

    :return: An array with three columns. The columns are
               0: latitude
               1: longitude
               2: response (unit following input R)
    """

    nx = len(ac)
    ny = len(al)
    W = np.zeros([nx * ny, 3])
    i = 0
    for ix in range(nx):
        for iy in range(ny):
            W[i, 0], W[i, 1] = latlon_one_antenna_point(
                lat0, lon0, ac[ix], al[iy], aa2sat, is_conical, re)
            W[i, 2] = Rg[ix, iy]
            i += 1

    return W


def inttime2angle(
    rpm: float,
    dt: float,
    is_conical: bool = False,
    za: float = 135
) -> float:
    """
    Calculates angular movement during an integration time.

    :param rpm:        Rotations per minute of the sensor
    :param dt:         Integration time
    :param is_conical: Cross-track (False, defualt) or concical (True) scanner.
    :param za:         Zenith angle. Only considered for conical scanners.
                       This is the angle at the sensor, not to be confused
                       with the incidence angle at the ground.

    :return: The angular movement.
    """

    dangle = (360 * rpm / 60) * dt
    if is_conical:
        dangle *= np.sin(np.deg2rad(za))

    return dangle


def latlon_dist_bearing(
    lat1: float,
    lon1: float,
    l: float,
    aa: float,
    re: float = 6378.1e3
    ): #-> tuple[float, float]:
    """
    Calculates position at distance and bearing.

    :param lat1: Latitude of start position.
    :param lon1: Longitude of start position.
    :param l:    Distance.
    :param aa:   Bearing, as azimuth angle.
    :param re:   Planet radius.

    :return: Tuple with lat and lon of new postion.
    """

    # Code from http://www.movable-type.co.uk/scripts/latlong.html
    # (but with short-cuts, such as asin(sin(lat2)) = lat2)
    dang = l / re
    cosdang = np.cos(dang)
    sindang = np.sin(dang)
    latrad = np.deg2rad(lat1)
    coslat = np.cos(latrad)
    sinlat = np.sin(latrad)
    aarad = np.deg2rad(aa)

    lat2 = sinlat * cosdang + coslat * sindang * np.cos(aarad)
    lon2 = lon1 + np.rad2deg(np.arctan2(np.sin(aarad) * sindang * coslat,
                                        cosdang - sinlat * lat2))
    lat2 = np.rad2deg(np.arcsin(lat2))

    return lat2, lon2


def latlon_one_antenna_point(
    lat0: float,
    lon0: float,
    ac: float,
    al: float,
    aa2sat: float,
    is_conical: bool = False,
    re: float = 6378.1e3
    ): #-> tuple[float, float]:
    """
    Calculates position for one antenne point given in (x,y).

    :param lat0:       Bore-sight latitude.
    :param lon0:       Bore-sight longitude.
    :param ac:         Across-track distance (at absolute aa=90).
    :param al:         Along-track distance (at absolute aa=90).
    :param aa2sat:     Azimuth angle, from ground to satellite [-180, 360].
    :param is_conical: Cross-track (False, default) or concical (True) scanner.
    :param re:         Planet radius.

    :return: Tuple with lat and lon of antenna point.
    """

    r = np.sqrt(ac * ac + al * al)
    aa = np.rad2deg(np.arctan2(ac, al))
    if is_conical:
        aa += aa2sat + 90
    else:
        if aa2sat <= 0 or aa2sat >= 180:
            aa += aa2sat + 90
        else:
            # Cross-track to left side
            aa += aa2sat - 90

    return latlon_dist_bearing(lat0, lon0, r, aa, re)


def length_from_ground(
    zp: float,
    ia: float,
    re: float = 6378.1e3,
    zg: float = 0.0
    ): #-> tuple[float, float]:
    """
    Find intersection with flat ground, for spherical planet.

    The function uses the incidence angle at the ground. If you instead have
    the zenith angle at the satellite, use instead length_to_ground.

    :param zp: Platform/satellite altitude.
    :param za: Incidence angle (at ground).
    :param re: Planet radius.
    :param zg: Ground altitude.

    :return: Tuple with distance to satellite and zenith angle at it.
    """

    if ia == 0:
        return zp - zg, 180
    else:
        beta = np.deg2rad(180-ia)
        q = (re + zp)/(np.sin(beta))
        na_rad = np.arcsin((re + zg)/q)
        return q * np.sin(beta + na_rad), 180-np.rad2deg(na_rad)


def length_to_ground(
    zp: float,
    za: float,
    re: float = 6378.1e3,
    zg: float = 0.0
    ): #-> tuple[float, float]:
    """
    Find intersection with flat ground, for spherical planet.

    The function uses the angle at the satellite. If you instead have
    the incidence angle at the ground, use instead length_to_ground.

    :param zp: Platform/satellite altitude.
    :param za: Zenith angle (at satellite).
    :param re: Planet radius.
    :param zg: Ground altitude.

    :return: Tuple with distance to ground and incidence angle.
    """

    if za == 180:
        return zp - zg, 0
    else:
        na_rad = np.deg2rad(180-za)
        q = (re + zg)/(np.sin(na_rad))
        ia = np.arcsin((re + zp)/q)
        return q * np.sin(ia - na_rad), np.rad2deg(ia)


def movement_add(
    x: np.ndarray,
    y: np.ndarray,
    R: np.ndarray,
    l: float,
    is_conical: bool = False,
    nsamples: int = 11
) -> np.ndarray:
    """
    Adds the effect of a movement.

    The movement is likely from scanning during integration time. That
    is, the function determined EFOV from IFOV. If l <= 0, nothing is done.

    Direction of scanning determined by scanning type. Conical
    broadens in the y-direction, and cross-track in x-direction.

    :param x:          x-grid (za or ac).
    :param y:          y-grid (aa or al).
    :param R:          Input response (IFOV).
    :param l:          Length of movement, in same unit as x and y.
    :param is_conical: Cross-track (False, defualt) or concical (True) scanner.
    :param nsamples:   Number of discrete positions inside movement to use.

    :return: Broadened response (EFOV).
    """

    # Nothing to do if l <= 0
    if l <= 0:
        return R

    # Shifts to apply
    dl = l / nsamples
    shifts = np.linspace((-l+dl)/2, (l-dl)/2, nsamples)

    # Variable for summed responses
    S = np.zeros_like(R)

    if is_conical:
        for i in range(len(x)):
            F = interpolate.interp1d(
                y,
                R[i, :],
                kind='cubic',
                fill_value="extrapolate"
            )
            for s in range(nsamples):
                S[i, :] += F(y + shifts[s])
    else:
        for i in range(len(y)):
            F = interpolate.interp1d(
                x,
                R[:, i],
                kind='cubic',
                fill_value="extrapolate"
            )
            for s in range(nsamples):
                S[:, i] += F(x + shifts[s])

    return S / nsamples


def norm_basic(
    x: np.ndarray,
    y: np.ndarray,
    R: np.ndarray
) -> np.ndarray:
    """
    Fast normalisation of integrated response.

    Performs a basic normalisation, assuming that both grids are
    equidistant. That is, the integrated response is normalised to
    one, assuming locally constant values.

    :param x: x-grid.
    :param y: y-grid.
    :param R: Input response.

    :return: Normalised response.
    """

    return R / (R.sum() * (x[1] - x[0]) * (y[1] - y[0]))


def norm_sum(
    R: np.ndarray
) -> np.ndarray:
    """
    Normalisation of summed values.

    :param R: Input response.

    :return: Normalised response.
    """

    return R / R.sum()


def scale2ground(
    za: np.ndarray,
    aa: np.ndarray,
    R:  np.ndarray,
    l:  float,
    ia: float): #-> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Scales footprint, from degress to distances at ground.

    :param za: Zenith angle grid
    :param aa: Azimuth angle grid
    :param R:  Input response.
    :param l:  Distance to ground.
    :param ia: Incidence angle.

    :return: Tuple ac grid, al grid and response
    """
    yf = (np.pi / 180) * l
    xf = yf / np.cos(np.deg2rad(ia))

    return xf*za, yf*aa, R/(xf*yf)

def get_efov(
    int_angle: float = 3.33
    ): #-> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    :return: Tuple with za, aa and R
    """
    #
    fwhm = 3.33
    x_fac = 3.5
    y_fac = 2
    dxy = 0.02
    #
    za, aa, R = gaussian2d(
        fwhm_x=fwhm,
        fwhm_y=fwhm,
        width_x=x_fac*fwhm,
        width_y=y_fac*fwhm,
        dx=dxy,
        dy=dxy)
    R = movement_add(
        x=za,
        y=aa,
        R=R,
        l=int_angle,
        is_conical=False)

    return za, aa, R


def get_land_sea_mask(
):
    """
    :return: Tuple with lat, lon, mask
    """
    #
    A = np.load("arome_arctic_lsm_2019100118.npz")
    #
    return {"lat": A['arr_1'].flatten(),
            "lon": A['arr_2'].flatten(),
            "value": A['arr_0'].flatten()}


def interp_lsmask(
        LSmask,
        lats,
        lons
):
    return griddata(points=(LSmask["lat"], LSmask["lon"]),
                    values=LSmask["value"],
                    xi=(lats, lons),
                    method='nearest')

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

# Antenna data
fwhm = 3.33
int_angle = 3.33
x_fac = 2.5
y_fac = 2
dxy = 0.05

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

#
# Calculations
#

# Create circular IFOV response
za, aa, R = gaussian2d(
    fwhm_x=fwhm,
    fwhm_y=fwhm,
    width_x=x_fac*fwhm,
    width_y=y_fac*fwhm,
    dx=dxy,
    dy=dxy)

# Go to EFOV - MM:Not now to compare first Antenna code and Footprint operator
R = movement_add(
    x=za,
    y=aa,
    R=R,
    l=int_angle,
    is_conical=False)

# Scale to ground
l, _ = length_from_ground(zsat, zenith)
ac, al, Rg = scale2ground(za, aa, R, l, zenith)

# Place on ground and plot
W = geopos_response(ac, al, Rg, lat, lon, azimuth, False)
W[:, 2] = norm_sum(W[:, 2])
#
#Compute footprint operator points on the ground
flat, flon = footprint(lat, lon, azimuth, scanp, zsat/1e3, sens)
#Less contours
#matplotlib.colormaps[name]`` or ``matplotlib.colormaps.get_cmap(obj)
cmap = plt.cm.get_cmap("Greys",10)
#plt.scatter(W[:, 1], W[:, 0], c=W[:, 2], s=5, cmap=cmap)
map.scatter(W[:, 1], W[:, 0], c=W[:, 2], latlon=True, s=50, cmap=cmap)
map.scatter(flon, flat, latlon=True, s=30, color='r', marker='.')
#plt.xlabel("Longitude [deg]")
#plt.ylabel("Latitude [deg]")
fwhm_x, fwhm_y = fwhm_by_interp(ac, al, Rg)
print(f"Antenna-code-Along-track: {fwhm_y/1e3:.2f} {'km'}")
print(f"Antenna-code-Across-track: {fwhm_x/1e3:.2f} {'km'}")
print(f"Antenna-code-Area: {((fwhm_x/1e3)/2.0)*((fwhm_y/1e3)/2.0)*np.pi:.2f} {'km2'}")

plt.title("Footprint, antenna pattern, and operator points on the ground \n Sensor: "+sensname+", Lat:"+str(lat)+", Lon:"+str(lon)+"")
plt.savefig(""+sensname+"_antenna_sampling.png")
