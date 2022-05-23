% This code is based on the function ilm_crystal_build_base.m, by I. Lobato
%
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


function base = tfm_crystal_build_base(uc, transformations)
    
    n_uc = size(uc, 1);
    base = zeros(256, size(uc, 2));
    
    ic = 0;
    ee2 = (1e-4)^2;

    for i_g = 1:size(transformations,1)
        [R,t] = tfm_transformation2Matrix(transformations(i_g,:));
        R = permute(R, [2, 1, 3]);
        for it = 1:n_uc
            r = uc(it, 2:4)*R + t;
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
    base = base(1:ic, :);
%     tfm_plot_crystal(base,'h',gca); view([0 0 1]);
end