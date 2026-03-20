# satellite-bufr-pp

A repository to alter or make certain processes with satellite BUFR data. It helps to obtain single observations, footprint representation (coming later), and else. The content(scripts, programs, else) is tailed for the use on ECMWF's HPC Atos.

## Single observation (directory: singleobs)

### Select an active radiance record/observation

Run your Harmonie-AROME assimilation and get a CCMA ODB with active radiance observations. Extract Obsvalue, latitude, longitude, channel number from the ODB and select the desired observation to be used as the single observation.
An example (*can be done differently*).

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
./driver.sh arg1 arg2 arg3 arg4 arg5 arg6
``` 

where

- arg1: Input satellite BUFR
- arg2: The name of the instrument (available: AMSU-A,MHS,ATMS,AWS)
- arg3: Observation value of the desired observation record in two decimal/digits
- arg4: Latitude of the desired observation record in radian
- arg5: Longitude of the desired observation record in radian
- arg6: Desired channel or channels separated by comma (e.g., 4 or 4,1,2). The first listed channel has to match the observation, lat, and lon values!

For example:

```bash
./driver.sh /path/to/your/scratch/hm_home/your_exp/archive/observations/year/month/day/hour/amsub.thinned.yourdatehour MHS 245.95 1.2266941927887 -0.0010541788682046 4 
```

or

```bash
./driver.sh /path/to/your/dir/awsyourdatehour AWS 215.32 1.1743936564234 0.42650261865135 7,4,5,12,14,15,16,17,18
```

#### Output files

The driver script create two BUFR files as output. A BUFR file which contains all channels for the given single obs as single subset. Another one which consists of only selected channels according to the list provided in arg6.

```bash
single_subset_$arg2.bufr
single_selected_$arg2.bufr
```

#### Check your output (optional, not manadatory)

```bash
module reset
module unload odb_api
module load python3 ecmwf-toolbox
bufr_dump -p $output_single_obs.bufr
```

#### Check the single obs on map

The driver script by default makes a plot about the location of the single obs
See the above example here:

![alt text](https://github.com/matemile/satellite-bufr-pp/blob/main/singleobs/mhs_example_location.png?raw=true)

and another with AWS radiance here:

![alt text](https://github.com/matemile/satellite-bufr-pp/blob/main/singleobs/aws_example_location.png?raw=true)

### Rerun your assimilation
Modify your Bator script and fetch the single obs BUFR

For example in scr/Bator

>
> #ln -sf $OBDIR/amsub.thinned$DTG ./BUFR.amsub
>
> ln -sf /your/git/satellite-bufr-pp/your-single-obs.bufr ./BUFR.amsub
>

## Footprint representation (directory: footprint)

After the single observation BUFR file is created, there's an opportunity to visualise the FOV and antenna pattern on the ground given by the sensor we have. So far this visualisation works with AMSU-A, MHS, and IASI radiances.

### Use the python program with proper arguments

```bash
python3 sampling_repres.py arg1 arg2 arg3
```

where

- arg1: Input satellite single obs BUFR (e.g., ../singleobs/single\_subset_$SENSOR.bufr)
- arg2: The sensorID of the IAL radiance coding (e.g., MHS = 15, IASI = 16, etc)
- arg3: The name of the sensor (e.g., MHS)

For example:

```bash
python3 sampling_repres.py ../singleobs/single_subset_MHS.bufr 15 MHS
```

### Check the generated plot

The python program produces a plot about the antenna pattern and possible footprint operator representation.

An exmaple for MHS single obs presented above:

![alt text](https://github.com/matemile/satellite-bufr-pp/blob/main/footprint/MHS_antenna_sampling.png?raw=true)
