#!/usr/bin/env python3
import os
import sys
import math
import numpy as np
from eccodes import (
    codes_bufr_new_from_file,
    codes_set,
    codes_set_array,
    codes_get_array,
    codes_get,
    codes_write,
    codes_release,
)

def isclose(a, b, atol):
    return abs(a - b) <= atol

def find_subsets(
    filename,
    target_lat=None,
    target_lon=None,
    lat_key="latitude",
    lon_key="longitude",
    extra_conditions=None,
    lat_tol=1e-2,
    lon_tol=1e-2,
):
    """
    Scan all BUFR messages and subsets in 'filename' to find subset indices
    that match the given conditions.

    extra_conditions: list of dicts like
        [{"key": "#4#brightnessTemperature", "value": 248.38, "tol": 1e-2},
    Keys can be array-valued or scalar; if array-valued, the same index as lat/lon is used.
    """
    if extra_conditions is None:
        extra_conditions = []

    results = []

    with open(filename, "rb") as fin:
        msg_index = 0
        while True:
            h = codes_bufr_new_from_file(fin)
            if h is None:
                break
            msg_index += 1

            # unpack compressed data into arrays
            codes_set(h, "unpack", 1)

            # latitude / longitude arrays (if present)
            try:
                lats = list(codes_get_array(h, lat_key))
                lons = list(codes_get_array(h, lon_key))
            except Exception:
                # no lat/lon in this message
                codes_release(h)
                continue

            nsub = len(lats)
            if len(lons) != nsub:
                # inconsistent, skip
                codes_release(h)
                continue

            # Prepare extra arrays or scalars
            extra_data = []
            for cond in extra_conditions:
                key = cond["key"]
                try:
                    arr = list(codes_get_array(h, key))
                    extra_data.append(("array", key, arr, cond["value"], cond.get("tol", 0.0)))
                except Exception:
                    # try scalar
                    try:
                        val = codes_get(h, key)
                        extra_data.append(("scalar", key, val, cond["value"], cond.get("tol", 0.0)))
                    except Exception:
                        # key not found; this message cannot match
                        extra_data = None
                        break

            if extra_data is None:
                codes_release(h)
                continue

            # Search subsets
            for i in range(nsub):
                ok = True

                #print('isclose',i,isclose(lats[i], target_dlat, lat_tol))
                if target_lat is not None and not isclose(lats[i], target_lat, lat_tol):
                    continue
                if target_lon is not None and not isclose(lons[i], target_lon, lon_tol):
                    continue

                for kind, key, data_val, target_val, tol in extra_data:
                    if kind == "array":
                        v = data_val[i]
                    else:
                        v = data_val  # scalar, same for all subsets in message
                    # numeric vs non-numeric
                    if isinstance(v, (int, float)) and isinstance(target_val, (int, float)):
                        if not isclose(float(v), float(target_val), tol):
                            ok = False
                            break
                    else:
                        if v != target_val:
                            ok = False
                            break

                if ok:
                    subset_index = i + 1  # 1-based
                    results.append(
                        {
                            "message": msg_index,
                            "subset": subset_index,
                            "lat": lats[i],
                            "lon": lons[i],
                        }
                    )

            codes_release(h)

    return results

def blank_unwanted_channels(infile, outfile, chlist, keep_channel):
    with open(infile, "rb") as fin, open(outfile, "wb") as fout:
        h = codes_bufr_new_from_file(fin)
        if h is None:
            raise RuntimeError("No message found")

        codes_set(h, "unpack", 1)
        
        for ch_num in chlist:
            if ch_num == keep_channel:
                print(f"-- Keeping channel {ch_num}")
                continue
            bt_key = f"#{ch_num}#brightnessTemperature"
            try:
                # Get length first
                bt_vals = codes_get_array(h, bt_key)
                n_vals = len(bt_vals)
                print(f"Original {bt_key}: {bt_vals[0]}")
                
                # Create blank values
                blank_vals = [1.0] * n_vals  # 1K = impossible BT
                codes_set_array(h, bt_key, blank_vals)
                print(f"Set {bt_key} to 1.0K")
            except Exception as e:
                print(f"Error {bt_key}: {e}")

        # CRITICAL: Pack the changes back into data section
        codes_set(h, "pack", 1)
        
        codes_write(h, fout)
        codes_release(h)
        print(f"Wrote to {outfile}")

if __name__ == "__main__":
    if len(sys.argv) < 8:
        print(f"Usage: {sys.argv[0]} BUFR_FILE VALUE LAT LON CH OUTFILE SENSORNAME", file=sys.stderr)
        print(
            f"{sys.argv[0]} BUFR_FILE VALUE LAT LON CH OUTFILE SENSORNAME",
            file=sys.stderr,
        )
        sys.exit(1)

    bufr_file = sys.argv[1]
    target_val = float(sys.argv[2])
    target_lat = float(sys.argv[3])
    target_lon = float(sys.argv[4])
    target_ch = int(sys.argv[5])
    out_file = sys.argv[6]
    sens = str(sys.argv[7])

    # Example: search specifically for your ATOVS channel 4 case
    extra_conditions = [
        {
            "key": "#"+str(target_ch)+"#brightnessTemperature",
            "value": target_val,
            "tol": 1e-2,
        },
    ]

    target_dlat = np.round(np.rad2deg(target_lat), 2)
    target_dlon = np.round(np.rad2deg(target_lon), 2)

    INSTRUMENT_CHANNELS = {
        'AMSU-A': list(range(1, 16)),        # 15 channels
        'MHS':    list(range(1, 6)),         # 5 channels
        'ATMS': list(range(1, 23)),        # 22 channels
        'AWS': list(range(1, 20)),        # 19 channels
        # Add more instruments as needed
    }

    if sens not in INSTRUMENT_CHANNELS:
        raise ValueError(f"Unknown instrument '{sens}'. Available: {list(INSTRUMENT_CHANNELS.keys())}")

    channels = INSTRUMENT_CHANNELS[sens]
    if target_ch not in channels:
        raise ValueError(f"Target channel {target_ch} not in {sens} channels {channels}")

    print(f"Instrument: {sens}, keeping channel {target_ch}, blanking {len(channels)-1} others")

    matches = find_subsets(
        bufr_file,
        target_lat=target_dlat,
        target_lon=target_dlon,
        extra_conditions=extra_conditions,
    )

    temp_file = "subset.bufr"

    with open(bufr_file, "rb") as fin, open(temp_file, "wb") as fout:
        msg_idx = 0
        while True:
            h = codes_bufr_new_from_file(fin)
            if h is None:
                break
            msg_idx += 1
            
            codes_set(h, "unpack", 1)
            
            if msg_idx == matches[0]["message"]:
                codes_set(h, "extractSubset", matches[0]["subset"])
                codes_set(h, "doExtractSubsets", 1)
                codes_write(h, fout)
                break
                
            codes_release(h)

    print(f"Extracted message {matches[0]['message']}, subset {matches[0]['subset']}")

    blank_unwanted_channels(temp_file, out_file, channels, keep_channel=target_ch)

    if os.path.isfile(temp_file):
        os.remove(temp_file)
    else:
        # If it fails, inform the user.
        print("Error: %s file not found" % temp_file)
