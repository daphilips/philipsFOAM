#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
makeType=${1:-libso}
set -x

echo $makeType
