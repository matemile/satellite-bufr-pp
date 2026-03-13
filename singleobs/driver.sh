#!/bin/bash

module reset
module unload odb_api
module load python3 ecmwf-toolbox

# Argument validation check
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <Input BUFR file> <Obsvalue in 2digits> <Latitude in radian>"
    echo "<Longitude in radian> <Channel number> <Output Singleobs BUFR file> <Instrument>"
    exit 1
fi
# Execution
echo "---------------";
echo "Input BUFR: $1";
echo "Obsvalue: $2";
echo "Latitude: $3";
echo "Longitude: $4";
echo "Channel num: $5";
echo "Output Singleobs BUFR: $6";
echo "Instrument: $7";
echo "---------------";

#python3 extract_write.py amsub.thinned2025092615 4 70.94 -7.65 248.38 single-out.bufr
# Extract and write chosen single obs in another BUFR
python3 extract_write.py $1 $2 $3 $4 $5 $6 $7

# Plot single obs on a map
python3 plot_singleo.py $6 $5
