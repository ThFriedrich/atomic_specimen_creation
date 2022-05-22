% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function [txyz_c, a_r, b_r, c_r, R] = tfm_align_z(atoms, a, b, c, alpha, beta, gamma, T)
    % Aligns the z-axis of the crystal with the vector T.
    % atoms = [atom_type, x, y, z]
    % a, b, c, alpha, beta, gamma = cell parameters
    % T = [x, y, z]
    % Returns:
    % txyz_c = [atom_type, x, y, z], aligned crystal
    % a_r, b_r, c_r = cartesian cell vectors of the aligned crystal
    % R = rotation matrix
    
    xyz = atoms(:,2:4);
    at_Z = atoms(:,1);
    
    g = tfm_direct_structure_matrix(a, b, c, alpha, beta, gamma);
    
    X = [1 0 0]; 
    Y = [0 1 0];
    Z = [0 0 1];
    
    A = X * g;

    T = T * g;
    T_rot = T;

    dr = fcn_rot_dir('y',T_rot);
    Txz = fcn_proj_XZ(T);
    ry = dr * fcn_vec_angle(Z,Txz);
    Rot_y = roty(ry);
    T_rot = T_rot * Rot_y;

    dr = fcn_rot_dir('x',T_rot);
    Tyz = fcn_proj_YZ(T_rot);
    rx = dr * fcn_vec_angle(Z,Tyz);
    Rot_x = rotx(rx);
    T_rot = T_rot * Rot_x;

    A_rot = A * Rot_y * Rot_x;
    dr = fcn_rot_dir('z',A_rot);
    Axy = fcn_proj_XY(A_rot);
    rz = dr * fcn_vec_angle(X,Axy);
    Rot_z = rotz(rz);
    T_rot = T_rot * Rot_z;

    R = g * Rot_y * Rot_x * Rot_z;
    xyz_R = xyz * R;
    
    txyz_c = [at_Z xyz_R];
    a_r = X * R;
    b_r = Y * R;
    c_r = Z * R;
end

function n = fcn_rot_dir(ax,vec)
    
    switch ax
        case 'x'
            vec(1) = [];
            vec(1) = vec(1) * -1;
        case 'y'
            vec(2) = [];
        case 'z'
            vec(3) = [];
    end
    b_vec = sign(vec);
    b_vec(b_vec==0) = 1;
    if isequal(b_vec,[1 1]) % Q1
        n = 1;
    elseif isequal(b_vec,[-1 1]) % Q2
        n = -1;  
    elseif isequal(b_vec,[1 -1]) % Q3
        n = -1;  
    elseif isequal(b_vec,[-1 -1]) % Q4
        n = 1;
    end
end

function theta = fcn_vec_angle(u,v)
    theta = atan2d(norm(cross(u,v)),dot(u,v));
end

function P = fcn_proj_YZ(Vec)
    P = Vec;
    P(1) = 0;
end

function P = fcn_proj_XZ(Vec)
    P = Vec;
    P(2) = 0;
end

function P = fcn_proj_XY(Vec)
    P = Vec;
    P(3) = 0;
end

function P = fcn_proj(Vec,Dir)
    P = dot(Dir,Vec)*Dir;
end
