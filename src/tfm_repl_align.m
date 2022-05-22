
function [atoms_collect, R] = tfm_repl_align(crystal_par, T_hkl, na, nb, nc, lx, ly, lz, b_plot, b_hkl, h)
    
    % Build crystal
    if isempty(crystal_par.transformations)
        atoms = ilm_crystal_build_base(crystal_par);
    else
        atoms = tfm_crystal_build_base(crystal_par.asym_uc, crystal_par.transformations);
    end
    
    % (HKL) -> [uvw]
    if b_hkl
        g_m = tfm_metric_tensor(a, b, c, alpha, beta, gamma);
        T_uvw = T_hkl / g_m;
    else
        T_uvw = T_hkl;
    end
    
    % Align vector T_uvw to cartesian Z-Axis
    [atoms, A, B, C, R] = tfm_align_z(atoms, a, b, c, alpha, beta, gamma, T_uvw);
    atoms_collect = atoms;

    % Duplicate along lattice vectors
    atoms_collect = tfm_loop_dim(atoms_collect, A, na);
    atoms_collect = tfm_loop_dim(atoms_collect, B, nb);
    atoms_collect = tfm_loop_dim(atoms_collect, C, nc);

    % Shift centre to zero , cut to box
    cm = mean(atoms_collect(:,2:4));
    atoms_collect(:,2:4) = atoms_collect(:,2:4) - cm;
    b_x = atoms_collect(:,2)<(-lx/2) | atoms_collect(:,2)>(lx/2); 
    b_y = atoms_collect(:,3)<(-ly/2) | atoms_collect(:,3)>(ly/2); 
    b_z = atoms_collect(:,4)<(-lz/2) | atoms_collect(:,4)>(lz/2); 
    atoms_collect(b_x|b_y|b_z,:) = [];
    [atoms_collect,sft] = ilm_spec_recenter(atoms_collect,lx,ly,lz);
    
    ref_xyz = max(-cm+sft,min(atoms_collect(:,2:4)));
    if b_plot
      ti = join([regexprep(formula,'\d+','_\{$0\}'), ' - [', join(string(T_hkl),' ') ']']);
      tfm_plot_crystal(atoms_collect, 'g', [R ref_xyz'], 'title', ti,'h', h)
    end
end