# satellite-bufr-pp

A repository to alter or make certain processes with satellite BUFR data. It helps to obtain single observations, footprint representation (coming later), and else. The content(scripts, programs, else) is tailed for the use on ECMWF's HPC Atos.

## Single observation

### Select an active radiance record/observation

Run your Harmonie-AROME assimilation and get a CCMA ODB with active radiance observations. Extract Obsvalue, latitude, longitude, channel number from the ODB and select the desired observation to be used as the single observation.
An example (*can be extracted differently*).

```bash
cd /path/to/your/scratch/hm_home/your_exp/yourdate_hour/odb_ccma
module use /perm/hlam/apps/modulefiles/lmod
module load odb_api
odbdup -i CCMA -o ccma -l CCMA
cd ccma
odbsql -q "select obsvalue,lat,lon,vertco_reference_1 from hdr,body"
```

Output

> ...
>
> 237.86 1.2562600703175 0.3584993550059 3
>
> 252.02 1.2562600703175 0.3584993550059 4
>
> 263.17 1.2562600703175 0.3584993550059 5
>
> 238.85 1.2016975872416 0.14394952071674 3
>
> 240.38 1.2266941927887 -0.0010541788682046 3
>
> 249.95 1.2266941927887 -0.0010541788682046 4
>
> 253.36 1.2266941927887 -0.0010541788682046 5
>
> 245.75 1.1098164740997 0.017175785168876 3
>
> 263.71 1.1098164740997 0.017175785168876 4
>
> 274 1.1098164740997 0.017175785168876 5
>
> ...

Choose the one which is suitable for your experiment.

> **249.95 1.2266941927887 -0.0010541788682046 4**

### Setup/Clone this repository

```bash
git clone git@github.com:matemile/satellite-bufr-pp.git
```

### Use driver bash script and give proper arguments for your need

Run the driver script as the following.

```bash
./driver.sh arg1 arg2 arg3 arg4 arg5 arg6 arg7
``` 

where

- arg1: Input satellite BUFR
- arg2: Observation value of the desired observation record in two decimal/digits
- arg3: Latitude of the desired observation record in radian
- arg4: Longitude of the desired observation record in radian
- arg5: Channel number (intiger) of the desired record
- arg6: Output satellite BUFR include only the desired single obs.
- arg7: The name of the instrument (available: AMSU-A,MHS,ATMS,AWS)

For example:

```bash
./driver.sh /path/to/your/scratch/hm_home/your_exp/archive/observations/year/month/day/hour/amsub.thinned.yourdatehour 245.95 1.2266941927887 -0.0010541788682046 4 your-single-obs.bufr MHS
```

#### Check your output (optional, not manadatory)
```bash
module reset
module unload odb_api
module load python3 ecmwf-toolbox
bufr_dump -E"python" $output_single_obs.bufr > $output_single_obs.py
```

### Rerun your assimilation
Modify your Bator script and fetch the single obs BUFR

For example in scr/Bator

>
> #ln -sf $OBDIR/amsub.thinned$DTG ./BUFR.amsub
>
> ln -sf /your/git/satellite-bufr-pp/your-single-obs.bufr ./BUFR.amsub
>

## Footprint representation
