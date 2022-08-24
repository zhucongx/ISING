#!/usr/bin/env bash
if [[ "$1" == [Rr] ]]; then
    mode=R
elif [[ "$1" == [Dd] ]]; then
    mode=D
else
  read -p "mode ([R]elease/[D]ebug): " mode
fi

if [[ $mode == [Rr] ]]; then
    mode=Release
    echo Release mode
elif [[ $mode == [Dd] ]]; then
    mode=Debug
    echo Debug mode
else
    echo Try agin...
    exit
fi

rm -rf build
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE="$mode" -S ..
fi
make -j 12
cd ..
