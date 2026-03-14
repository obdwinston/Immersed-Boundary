#!/bin/bash
set -e

python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip -q
pip install -r requirements.txt -q

mkdir -p build outputs/data
rm -f outputs/data/fields_*.dat outputs/data/forces.dat

gfortran -O3 -march=native -fopenmp \
    -J build -I build -o main \
    src/config.f90 \
    src/solver.f90 \
    src/io.f90 \
    src/main.f90

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
export OMP_STACKSIZE=512M

./main

python animate.py
