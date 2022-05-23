% Copyright (c) 2010, Francisco de Castro
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of  nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% Source: scatter3sph (https://www.mathworks.com/matlabcentral/fileexchange/27112-scatter3sph), MATLAB Central File Exchange. Retrieved May 18, 2022. 
%
% Simplified and adapted the original version here

function h_pl = scatter3sph(h, XYZ, sp_size, color, varargin)

if nargin < 2
    sp_size = 0.1*max(XYZ);
end

if nargin < 3
    color = [0 0 1];
end

%-- Sphere facets
[sx,sy,sz]= sphere(32);

%-- Plot spheres
hold(h, 'on')
for j= 1:size(XYZ,1)
	h_pl = surf(h, sx*sp_size+XYZ(j,1), sy*sp_size+XYZ(j,2), sz*sp_size+XYZ(j,3),...
		'LineStyle','none',...
		'FaceColor',color,...
		'FaceAlpha',1);
end
view(h, [1 1 1])
light(h);
lighting(h,'gouraud');
hold(h, 'off')

