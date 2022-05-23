function crystal_par = tfm_get_uc_from_cif(cif_path)  
    
    [crystal_par.asym_uc, crystal_par.a, crystal_par.b, crystal_par.c,...
    crystal_par.alpha, crystal_par.beta, crystal_par.gamma,... 
    crystal_par.sgn, crystal_par.hmg, crystal_par.transformations,...
    crystal_par.formula] = tfm_import_cif(cif_path);
    
    if isempty(crystal_par.transformations)
        crystal_par.atoms = ilm_crystal_build_base(crystal_par);
    else
        crystal_par.atoms = tfm_crystal_build_base(crystal_par.asym_uc, crystal_par.transformations);
    end
end