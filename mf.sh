#!/bin/bash

# example usage: ./mf.sh ~/data/matchings/out.bibsonomy-2ti ~/data/matchings/results/bibsonomy > /dev/null
if [ "$#" -ne 2 ]; then
    echo "Wrong number of parameters"
    echo "Usage: ./mf.sh input_file output_dir"
    exit
fi
echo "input file = '$1'" 
echo "output dir = '$2'"
mkdir -p $2

if [ ! -f "$1" ]; then
    echo "Error: input file does not exist"
    exit
fi

if [ ! -f "./mf" ]; then
    g++ -O3 -DNDEBUG -DFULL_OUTPUT mf.cpp -o mf
fi
if [ ! -f "./bk" ]; then
    g++ -O3 -DNDEBUG bk.cpp -o bk
fi

if [ ! -f "$2/decomp" ]; then
    echo "$time ./mf < $1"
    ( time ./mf < $1 ) 2>&1 | tee -a "$2/log.aux"
    grep -v '#' "$2/log.aux" > "$2/decomp.log"
    grep '#' "$2/log.aux" | cut -d' ' -f2- > "$2/decomp"
    rm "$2/log.aux"
    rm "$2/matchings" 2> /dev/null
else    
    echo "decomp file exists; skipping..."
fi    

if [ ! -f "$2/matchings" ]; then
    echo "time ./bk $2/decomp" 
    ( time ./bk "$2/decomp" ) 2>&1 | tee -a "$2/matchings_log.aux"
    grep -v '#' "$2/matchings_log.aux" > "$2/matchings.log"
    grep '#' "$2/matchings_log.aux" | cut -d' ' -f2- > "$2/matchings"
    rm "$2/matchings_log.aux"
else
    echo "matchings file exists; skipping..."
fi
