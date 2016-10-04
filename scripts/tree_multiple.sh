#!/bin/bash

# Copyright 2016 Richard Rodrigues, Nyle Rodgers, Mark Williams, Virginia Tech
#
# This file is part of Coremic.
#
# Coremic is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Coremic is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Coremic. If not, see <http://www.gnu.org/licenses/>.
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
