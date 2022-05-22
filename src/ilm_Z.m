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

function[Z] = ilm_Z(str_Z)
    a_str_Z = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', ...
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', ...
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',...
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',...
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'};

    if(ischar(str_Z))
        str_Z = {str_Z};
    end
    str_Z = strtrim(str_Z);
    a_str_Z = strtrim(a_str_Z);
    n_str_Z = numel(str_Z);
    Z = zeros(n_str_Z, 1);
    for ik=1:n_str_Z
        Z(ik) = find(strcmpi(str_Z{ik}, a_str_Z));
    end
end