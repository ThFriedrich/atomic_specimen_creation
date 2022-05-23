% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function [atoms, R, crystal_par] = tfm_align_duplicate_cut(cif_path, T_hkl, rot_z, lx, ly, lz, b_plot, b_hkl, h)
    % path = path to cif file or crystal parameters loaded already
    % T_hkl = zone axis to align to cartesian z
    % na, nb, nc = number of unit cell replications along lattice vectors
    % lx, ly, lz = box size
    % b_plot = boolean, plot the atoms?
    % b_hkl = boolean, zone axis given in hkl, otherwise uvw?
    % h = figure handle

    if isstruct(cif_path)
        crystal_par = cif_path;
    else
        [crystal_par.asym_uc, crystal_par.a, crystal_par.b, crystal_par.c,...
        crystal_par.alpha, crystal_par.beta, crystal_par.gamma,... 
        crystal_par.sgn, crystal_par.hmg, crystal_par.transformations,...
        crystal_par.formula] = tfm_import_cif(cif_path);
    end
    

    if isempty(crystal_par.transformations)
        atoms = ilm_crystal_build_base(crystal_par);
    else
        atoms = tfm_crystal_build_base(crystal_par.asym_uc, crystal_par.transformations);
    end

    % (HKL) -> [uvw]
    if b_hkl
        g_m = tfm_metric_tensor(crystal_par.a, crystal_par.b, crystal_par.c,... 
            crystal_par.alpha, crystal_par.beta, crystal_par.gamma);
        T_uvw = T_hkl / g_m;
    else
        T_uvw = T_hkl;
    end
    
    % Align vector T_uvw to cartesian Z-Axis
    [atoms, A, B, C, R] = tfm_align_z(atoms, crystal_par.a, crystal_par.b, crystal_par.c,...
        crystal_par.alpha, crystal_par.beta, crystal_par.gamma, T_uvw, rot_z);

    % Duplicate along lattice vectors
    if all([lx,ly,lz]>0)
        [atoms, cm] = tfm_loop_dim(atoms, [A' B' C'], lx, ly, lz, false);
    else
        cm = [0 0 0];
    end

    at_rng = max(atoms(:,2:4))-min(atoms(:,2:4));
    [atoms,sft] = ilm_spec_recenter(atoms,at_rng(1),at_rng(2),at_rng(3));
    
    if b_plot
      ref_xyz = min(max(-cm+sft,min(atoms(:,2:4))),max(atoms(:,2:4)));
      ti = join([regexprep(crystal_par.formula,'\d+','_\{$0\}'), ' - [', join(string(T_hkl),' ') ']']);
      tfm_plot_crystal(atoms, 'g', [R ref_xyz'], 'title', ti,'h', h)
    end
end