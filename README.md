# satellite-bufr-pp
A repository to pre-process satellite BUFR data and to obtain single observations, footprint geometry(coming later), etc. The script and python programs can be used on ECMWF's HPC Atos. For the time being, it has been tested with MHS BUFR, but it is intended to make it generic and applicable for other satellite radiances as well.

# Select an active radiance record/observation
Run your Harmonie-AROME assimilation and get a CCMA ODB with active radiance observations. Extract Obsvalue, latitude, longitude, channel number from the ODB and select the desired observation to be used as the single observation.

# Setup/Clone this repository
```bash
git clone git@github.com:matemile/satellite-bufr-pp.git
```

# Use driver bash script and give proper arguments for your need
Run the driver script as the following.
```bash
./driver.sh arg1 arg2 arg3 arg4 arg5 arg6
``` 
where 
arg1 - Input satellite BUFR
arg2 - Channel number (intiger) of the desired record
arg3 - Latitude of the desired observation record in two decimal/digits
arg4 - Longitude of the desired observation record in two decimal/digits
arg5 - Observation value of the desired observation record in two decimal/digits
arg6 - Output satellite BUFR include only the desired single obs.

# Check your output (optional, not manadatory)
```bash
module reset
module unload odb_api
module load python3 ecmwf-toolbox
bufr_dump -E"python" $output_single_obs.bufr > $output_single_obs.py
```

# Rerun your assimilation
Modify Bator script and fetch the single obs output BUFR by Bator
