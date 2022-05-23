% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function g = tfm_metric_tensor(a, b, c, alpha, beta, gamma)

ab_cos_g = a * b * cosd(gamma);
ac_cos_b = a * c * cosd(beta);
bc_cos_a = b * c * cosd(alpha);

g = [a^2,       ab_cos_g    ac_cos_b;...
     ab_cos_g   b^2,     	bc_cos_a;...
     ac_cos_b   bc_cos_a 	c^2];
 