#!/bin/sh
F90=ifort
$F90 station.f90
./a.out > station_reg2.tbl
rm a.out
