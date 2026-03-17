import sys
import numpy as np
import matplotlib.pyplot as plt
from extract_write import parse_channels_arg
from mpl_toolkits.basemap import Basemap
from eccodes import (
    codes_bufr_new_from_file, codes_set, codes_get_array,
    codes_get, codes_release
)

def plot_single_obs_auto(bufr_file, channel_num):
    """
    Plot single obs on a Basemap with an automatically chosen domain
    around the observation.
    """
    # --- Read BUFR ---
    with open(bufr_file, "rb") as fin:
        h = codes_bufr_new_from_file(fin)
        if h is None:
            raise RuntimeError("No message in BUFR file")
        codes_set(h, "unpack", 1)

        lat = float(codes_get_array(h, "latitude")[0])
        lon = float(codes_get_array(h, "longitude")[0])
        bt_key = f"#{channel_num}#brightnessTemperature"
        bt_value = float(codes_get_array(h, bt_key)[0])

        year = codes_get(h, "year")
        month = codes_get(h, "month")
        day = codes_get(h, "day")
        hour = codes_get(h, "hour")
        minute = codes_get(h, "minute")
        satellite = codes_get(h, "satelliteIdentifier")

        codes_release(h)

    # choose projection – cylindrical is simple and robust
    fig = plt.figure(figsize=(8, 6))
    m = Basemap(width=1000000,height=1000000, resolution='i',projection="stere",
            lat_0=lat,lon_0=lon
    )

    # --- Draw map ---
    m.drawcoastlines(linewidth=1.0)
    m.drawparallels(np.arange(-90, 91, 5), labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180, 181, 10), labels=[0,0,0,1])

    x, y = m(lon, lat)
    m.scatter(x, y, c=bt_value, s=40, linewidth=0)

    plt.title(
        f"Single BUFR Observation\n"
        f"{year}-{month:02d}-{day:02d} {hour:02d}:{minute:02d} UTC | "
        f"SatID={satellite} | BT={bt_value:.2f}K | Ch{channel_num}",
        fontsize=11,
    )
    plt.tight_layout()
    fig.savefig('single_obs_location.png', dpi=300, format='png')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} SENSORNAME CHANNEL", file=sys.stderr)
        print(
            f"{sys.argv[0]} SENSORNAME CHANNEL",
            file=sys.stderr,
        )
        sys.exit(1)

    bufr_file = "single_selected_"+sys.argv[1]+".bufr" #sys.argv[1]
    allchannel = str(sys.argv[2])

    channel = parse_channels_arg(allchannel)

    print('bufr_file,channel',bufr_file,channel[0])

    # Usage
    plot_single_obs_auto(bufr_file, channel[0])
