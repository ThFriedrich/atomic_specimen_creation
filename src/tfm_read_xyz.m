% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function [atoms, lx, ly, lz] = tfm_read_xyz(varargin)
% Read atomic positions from xyz file and add rms3d to atom types
% First Input argument is the path to the xyz file.
% The rest are pairs of Atom description, rms3d value.
% e.g. [atoms, lx, ly, lz] = tfm_read_ap_xyz('/FeC_example.xyz','Fe',0.087,'C',0.091)
    
    % Open file
    path = varargin{1};
    fid = fopen(path, 'r');
    nh = 0;
    % Determine Header Lines
    while true
        str = fgetl(fid);
        C = textscan(str, '%s %f %f %f');
        if ~any(isnan(C{2})) && ~any(isempty(C{2}))
            break;
        end
        nh = nh + 1;  
    end
    % Read Data and close
    C = textscan(fid, '%s %f %f %f','CollectOutput',true,'HeaderLines',nh);
    fclose(fid);
    
    % Convert Label to Z and create rms3d Vector based on atom labels
    Z = ilm_Z(C{1});
    rms3d = zeros(size(Z,1),1);
    
     if(nargin<2)
        rms3d(:) = 0.085;
     else
        for id = 2:2:(nargin-1)
           b = char(C{1}) == char(varargin{id});
           rms3d(b) = varargin{id+1};
        end
     end

    xyz = C{2};
    xyz = xyz - min(xyz);
    atoms = [Z, xyz, rms3d, ones(size(Z,1),1)];
    [lx, ly, lz] = ilm_vect_assign(max(xyz));

end