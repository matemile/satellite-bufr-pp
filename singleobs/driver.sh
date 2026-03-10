#!/bin/bash

module reset
module unload odb_api
module load python3 ecmwf-toolbox

# Argument validation check
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <Input BUFR file> <Channel number> <Latitude in DEG 2digits>"
    echo "<Longitude in DEG 2digits> <Obsvalue in 2digits> <Output Singleobs BUFR file>"
    exit 1
fi
# Execution
echo "Input BUFR: $1";
echo "Channel Num.: $2";
echo "Latitude: $3";
echo "Longitude: $4";
echo "Obsvalue: $5";
echo "Output Singleobs BUFR: $6";

#python3 extract_write.py amsub.thinned2025092615 4 70.94 -7.65 248.38 single-out.bufr
python3 extract_write.py $1 $2 $3 $4 $5 $6
