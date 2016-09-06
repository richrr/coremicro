#!/bin/bash
INTEREST=$1
OUT=$2
OUTPUT=$3

# Directory of script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VISUALIZE_CORE="python2 $DIR/visualize_core.py"

RASTER="png"
VECTOR="svg"

$VISUALIZE_CORE $INTEREST $OUT $OUTPUT.$RASTER $OUTPUT.$VECTOR

for frac in 75 80 85 90 95 100; do
    $VISUALIZE_CORE -n $frac -x $frac $INTEREST $OUT ${OUTPUT}_$frac.$RASTER ${OUTPUT}_$frac.$VECTOR
done
