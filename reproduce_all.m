close all; clear; clc;

addpath("src/")
saveresults = true;

% generate 9 pictures (Fig 1 right, Fig. 5 - 7)
biomol_plots(saveresults)

% generate 2 pictures (Fig. 8) and 1 table (Tab. 4)
biomol_scales(saveresults)

% generate 2 tables (Tab. 3)
biomol_tabs()

% generate 5 pictures (Fig 1 left, Fig. 2 - 3)
newton_plots(saveresults)

% generate 2 tables (Tab. 1 - 2)
newton_tabs()

% generate 1 picture (Fig. 9)
lattice_plots(saveresults)
