% This file is part of Multem.
% Copyright 2022 Ivan Lobato <Ivanlh20@gmail.com>% 
% Multem is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version of the License, or
% (at your option) any later version.% 
% Multem is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.% 
% You should have received a copy of the GNU General Public License
% along with Multem. If not, see <http:// www.gnu.org/licenses/>.

function [atoms, sft] = ilm_spec_recenter(atoms, lx, ly, lz, opt)
    if(nargin<3)
        ly = lx;
    end
    
    if(nargin==3)
        lz = lx;
    end    
    
    if(nargin==3)
        opt = 2;
    elseif(nargin<5)
        opt = 1;
    end
    
    xyz_min = min(atoms(:, 2:4));
    atoms(:, 2:4) = atoms(:, 2:4)- xyz_min;
    xys_m = 0.5*([lx, ly, lz]-max(atoms(:, 2:4)));
    if(opt==1)
        atoms(:, 2:4) = atoms(:, 2:4) + xys_m;
        if(nargout==2)
            sft = -xyz_min + xys_m;
        end
    else
        atoms(:, 2:3) = atoms(:, 2:3) + xys_m(1:2);
        if(nargout==2)
            sft = -xyz_min(1:2) + xys_m(1:2);
        end
    end  
end