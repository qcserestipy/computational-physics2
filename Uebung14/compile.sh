#!/bin/bash

gfortran cp_mod.f95 FT_mod.f95 main.f95 && ./a.out && echo "succes!" && sleep 1
gnuplot nice.plot
evince out.pdf
