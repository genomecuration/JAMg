#!/usr/bin/env bash

set -e -v

cd htslib
mkdir -p build
autoheader
autoconf
./configure --prefix=`pwd`/build/
make
make install


