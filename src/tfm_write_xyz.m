% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function tfm_write_xyz(atoms,filename)
    % Write XYZ file
    %
    % INPUT
    %   atoms = [atom_type, x, y, z]
    %   filename: filename
    
    la = tfm_Z_str(atoms(:,1));
    xyz = atoms(:,2:4);
    n_at = size(xyz,1);
    fid = fopen (filename, "w");
    fprintf(fid, "%d\n",n_at);
    for ix = 1:n_at
        fprintf(fid, "%s\t%.4f\t%.4f\t%.4f\n", la{ix}, xyz(ix,:));
    end
    fclose (fid);
end
