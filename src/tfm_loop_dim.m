% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function atoms_collect = tfm_loop_dim(atoms_collect, V, n)
    atoms_init = atoms_collect;
    P = [atoms_init(:,2:end)'; ones(1,size(atoms_init,1))];    
    for ic = 1:n
        Tr = eye(4,4);
        Tr(1:3,4) = ic*V;
        
        tr_xyz = Tr*P;
        tr_xyz = tr_xyz(1:end-1,:)';
        
        atoms_dup = [atoms_init(:,1) tr_xyz];
        atoms_collect = cat(1, atoms_collect, atoms_dup);
    end
end