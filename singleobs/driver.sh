#!/bin/bash

module reset
module unload odb_api
module load python3 ecmwf-toolbox

# Argument validation check
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <Input BUFR file> <Instrument> <Obsvalue in 2digits> <Latitude in radian>"
    echo "<Longitude in radian> <Channel list> "
    exit 1
fi
# Execution
echo "---------------";
echo "Input BUFR: $1";
echo "Instrument: $2";
echo "Obsvalue: $3";
echo "Latitude: $4";
echo "Longitude: $5";
echo "Channel list: $6";
echo "---------------";

#python3 extract_write.py amsub.thinned2025092615 4 70.94 -7.65 248.38 single-out.bufr
# Extract and write chosen single obs in another BUFR
python3 extract_write.py $1 $2 $3 $4 $5 $6

# Plot single obs on a map
python3 plot_singleo.py $2 $6
