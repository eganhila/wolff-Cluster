#!/usr/local/bin/MathematicaScript -script

spins = Import["spins.dat"];
Export["spinPlot.png",MatrixPlot[spins]]
