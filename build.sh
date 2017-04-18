#!/bin/bash

module unload matlab
module load matlab/R2012b

mkdir ../bin

mcc -N -m -R '-nodisplay' -v split_tissue_maps.m -d ../bin -a ./NIfTI_20140122 -a ./TVAL3_v1.0 -a ./TVAL3D

