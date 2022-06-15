clear; clc;

% file_path = path to .cif file
file_path = 'VO_P42.cif';

%hkl Zone axis you want to see
% hkl = [3 1 0];
uvm = [1 1 1];

% lx, ly, lz = box size
lx = 40; ly = 40; lz = 10;

% Wanna see the structure plotted?
b_plot = true;

atoms = tfm_align_duplicate_cut(file_path, uvm, 0, lx, ly, lz, b_plot, false, gca);

% For the Multem input one has to add columns for rms3d and region
rms3d = ones(size(atoms,1),1) * 0.085;
atoms = [atoms, rms3d, ones(size(atoms,1),1)];



