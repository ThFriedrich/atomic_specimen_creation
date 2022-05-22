% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function g = tfm_direct_structure_matrix(a, b, c, alpha, beta, gamma)
  % direct structure matrix
  % transforms crystal coordinates to Cartesian coordinates.
  % M. De Graef, Introduction to Conventional Transmission Electron Microscopy, vol. 51, no. 2–3. 2003. p. 57
  
    cos_a = cosd(alpha);
    cos_b = cosd(beta);
    cos_g = cosd(gamma);    
    sin_g = sind(gamma);
    
    cos_bg = cos_b * cos_g;
    cos_bga = cos_bg - cos_a;
    
    v2 = a^2 * b^2 * c^2 * (1-cos_a^2-cos_b^2-cos_g^2+2*cos_a*cos_bg);
    v = sqrt(v2);
    
    g = [a      b*cos_g     c*cos_b;...
        0       b*sin_g    -c*cos_bga/sin_g;...
        0       0           v/(a*b*sin_g)];
    
    g = g.';
end