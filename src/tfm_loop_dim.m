% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function [atoms, cm_o] = tfm_loop_dim(atoms, V, lx, ly, lz, b_test)
    lxyz = sqrt(sum([lx ly lz].^2));
    for nd = 1:length(V)
        n_at = size(atoms,1);
        Tr = eye(4, 4,'single');
        ic = 1;
        while true       
            Tr(1:3,4) = ic*V(:,nd);
            tr_xyz = (Tr*[atoms(1:n_at,2:end)'; ones(1,n_at,'single')])';
            tr_xyz(:,end) = [];
            cm = fcn_center(atoms);
            atoms_it = [atoms(1:n_at,1) tr_xyz];
            atoms = cat(1, atoms, atoms_it);
            ic = ic + 1;
            if atoms_in(atoms_it, cm, lxyz) == 0
                break;
            end 
            if b_test
                tfm_plot_crystal([atoms(:,1) atoms(:,2:4)]); axis on; hold on;
                plotcube([lx, ly, lz],cm-[lx, ly, lz]./2,0.5,'y'); hold off;
            end
        end
    end
    
    % Shift centre to zero , cut to box
    [cm, cm_o] = fcn_center(atoms);
    atoms(:,2:4) = atoms(:,2:4) - cm;
    [~, b_out] = atoms_crop(atoms, [lx ly lz]);
    
    if b_test
        plot_in_out(atoms,b_out);
    end   
    
    atoms(b_out,:) = [];
    
end

function [cm, cm_o] = fcn_center(atoms)
    cm =  mean(atoms(:,2:4),1);
    cm_o = cm;
    [~, z0_i] = min(abs(atoms(:,4)-cm(3)));
    cm(3) = cm(3) - (atoms(z0_i,4)-cm(3));
end

function [at_in, out] = atoms_in(atoms, cm, lxyz)
    l = lxyz/2;
    out =   atoms(:,2)<(-l+cm(1)) | atoms(:,2)>l+cm(1) |... 
            atoms(:,3)<(-l+cm(2)) | atoms(:,3)>l+cm(2) |...  
            atoms(:,4)<(-l+cm(3)) | atoms(:,4)>l+cm(3); 
    at_in = sum(~out);
end

function [at_in, out] = atoms_crop(atoms, lxyz)
    l = lxyz/2;
    out =   atoms(:,2)<(-l(1)) | atoms(:,2)>l(1) |... 
            atoms(:,3)<(-l(2)) | atoms(:,3)>l(2) |...  
            atoms(:,4)<(-l(3)) | atoms(:,4)>l(3); 
    at_in = sum(~out);
end

function plot_in_out(atoms,b_dist)
    figure(6); clf;
    subplot(1,3,1)
    scatter(atoms(:,2),atoms(:,3),'xr'); hold on;
    scatter(atoms(~b_dist,2),atoms(~b_dist,3),'ob'); 
    scatter(0,0,'k*'); hold off;
    axis equal; title('x-y');

    subplot(1,3,2)
    scatter(atoms(:,2),atoms(:,4),'xr'); hold on;
    scatter(atoms(~b_dist,2),atoms(~b_dist,4),'ob'); 
    scatter(0,0,'k*'); hold off;
    axis equal; title('x-z');

    subplot(1,3,3)
    scatter(atoms(:,3),atoms(:,4),'xr'); hold on;
    scatter(atoms(~b_dist,3),atoms(~b_dist,4),'ob'); 
    scatter(0,0,'k*'); hold off;
    axis equal; title('y-z');
end