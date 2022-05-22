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
%
% This function needs the binary file 'space_groups.mat' which was also authored by
% Ivan Lobato.

function[base]=ilm_crystal_build_base(crystal_par)

    sg = importdata('space_groups.mat');
    sym = sg(crystal_par.sgn);

    uc = crystal_par.asym_uc;
    n_uc = size(uc, 1);
    M = permute(sym.M, [2, 1, 3]); % tranpose fro x*M
    n_M = size(M, 3);
    n_tr_g = size(sym.tr_g, 1);
    base = zeros(256, size(uc, 2));
    ic = 0;
    ee2 = (1e-4)^2;
    cnt = 0;
    for i_g = 1:n_tr_g
        tr_g = sym.tr_g(i_g, :);
        for is = 1:n_M
            cnt = cnt + 1;
            M_ik = M(:, :, is);
            tr = sym.tr(is, :) + tr_g;
            tr = tr - floor(tr);
            for it = 1:n_uc
                r = uc(it, 2:4)*M_ik + tr;
                r = r - floor(r);
                if((ic==0))
                    ic = ic + 1;
                    base(ic, :) = [uc(it, 1), r, uc(it, 5:end)];
                else
                    ii = find(sum((base(1:ic, 2:4)-r).^2, 2)<ee2);
                    if(size(ii, 1)==0)
                        ic = ic + 1;
                        base(ic, :) = [uc(it, 1), r, uc(it, 5:end)];
                    end
                end
            end
        end
    end
    base = base(1:ic, :);
    figure(4); clf;
    tfm_plot_crystal(base,'h',gca); view([0 0 1]);
end