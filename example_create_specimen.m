clear; clc;

% file_path = path to .cif file
file_path = 'VO_P42.cif';

%hkl Zone axis you want to see
% hkl = [3 1 0];
uvm = [0 1 1];

% na, nb, nc = number of unit cell replications along lattice vectors
na = 20; nb = 20; nc = 20; 

% lx, ly, lz = box size
lx = 40; ly = 40; lz = 20;

% Wanna see the structure plotted?
b_plot = true;

atoms = tfm_align_duplicate_cut(file_path, uvm, na, nb, nc, lx, ly, lz, b_plot, false, gca);

% For the Multem input one has to add columns for rms3d and region
rms3d = ones(size(atoms,1),1) * 0.085;
atoms = [atoms, rms3d, ones(size(atoms,1),1)];



