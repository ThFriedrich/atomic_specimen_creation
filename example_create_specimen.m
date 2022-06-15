clear; clc;

% file_path = path to .cif file
file_path = 'VO_P42.cif';

%hkl Zone axis you want to see
% hkl = [3 1 0]; b_hkl = true;  % plane normal
uvm = [1 1 1]; b_hkl = false;   % lattice vector

% rotation angle between the cartesian x axis and the projected
% crystallographic vector a around the cartesian z axis
rot_z = 0;

% lx, ly, lz = box size
lx = 40; ly = 40; lz = 10;

% Wanna see the structure plotted?
b_plot = true;

atoms = tfm_align_duplicate_cut(file_path, uvm, rot_z, lx, ly, lz, b_plot, b_hkl);

% For the Multem input one has to add columns for rms3d and region
rms3d = ones(size(atoms,1),1) * 0.085;
atoms = [atoms, rms3d, ones(size(atoms,1),1)];



