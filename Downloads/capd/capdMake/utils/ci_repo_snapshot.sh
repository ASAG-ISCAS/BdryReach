#!/bin/bash

set -e


#it is quite difficult to set PATH for jenkins-slave started by launchd on OSX.
export PATH=/usr/local/bin:$PATH
source "$(dirname $0)/ci_configure_flags.sh"
source "$(dirname $0)/ci_funcs.sh"


output=$PWD/output
cd workdir/src

echo "Calling autoreconf"
autoreconf -fi -v --no-recursive # only in curent directory to generate capdMake/capd_version_number.raw
version=$(cat capdMake/capd_version_number.raw)

echo "Version is: $version"

mkdir -p "$output"
tar --transform "s/^\./capd-$version/" -czf $output/capd-${version}.tar.gz .

echo $(date +%Y%m%d_%H%M) > $output/build_date

echo "Archive created:"
ls -ltr $output
