#!/bin/bash

# Copyright 2016, 2017 Richard Rodrigues, Nyle Rodgers, Mark Williams,
# Virginia Tech
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
TABLE=$1
GROUPFILE=$2
OUTPUT_DIRECTORY=$3

# Directory of script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESIZE_CMD="python2 $DIR/resize.py"

SIZES=(2 10 25 50)

echo 1
$RESIZE_CMD $TABLE $GROUPFILE $OUTPUT_DIRECTORY/1x1.biom $OUTPUT_DIRECTORY/1x1.tsv -c 1 -r 1
for SIZE in ${SIZES[*]}; do
    echo $SIZE
    $RESIZE_CMD $TABLE $GROUPFILE $OUTPUT_DIRECTORY/${SIZE}x1.biom $OUTPUT_DIRECTORY/${SIZE}x1.tsv -c $SIZE
    $RESIZE_CMD $TABLE $GROUPFILE $OUTPUT_DIRECTORY/1x${SIZE}.biom $OUTPUT_DIRECTORY/1x${SIZE}.tsv -r $SIZE
    $RESIZE_CMD $TABLE $GROUPFILE $OUTPUT_DIRECTORY/${SIZE}x${SIZE}.biom $OUTPUT_DIRECTORY/${SIZE}x${SIZE}.tsv -c $SIZE -r $SIZE
done
