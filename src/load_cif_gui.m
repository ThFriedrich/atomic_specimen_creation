% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

% Create a new figure with given size and minimum size
close all;
screensize = get(0, 'Screensize');
s = round(screensize(3:4)*0.8);
hpf.fig = figure('units','pixels','outerposition',[screensize(3)/2-s(1)/2 screensize(4)/2-s(2)/2 s(1) s(2)],'Name','Cif Import','NumberTitle','off','Visible','on','Resize','on','MenuBar','none','ToolBar','figure');
figRSfun = @(~,~) set(hpf.fig, 'position', max([0 0 900 550], hpf.fig.Position));
hpf.fig.SizeChangedFcn = figRSfun;

% Normalized vertical panel split position and border width
v_sec = 0.15;
h_sec = 0.8;
br = 0.01;

% Create 3 sub-panels
% Image Panel & Axis
hpf.image.pan = uipanel('Parent',hpf.fig,'Title','Visualizations','units','normalized','Position',[br v_sec+br h_sec-br 1-v_sec-(3*br)],'ShadowColor',[0 0 0],'ForegroundColor',[0 0 0],'HighlightColor',[0.95 0.95 0.95],'BackgroundColor',[0.8 0.8 0.8]);
hpf.image.ax1 = subplot(1,2,1,'Parent',hpf.image.pan);
hpf.image.ax2 = subplot(1,2,2,'Parent',hpf.image.pan);

v_sec2 = v_sec*4;
% Atomic Info Panel
hpf.inf.panA = uipanel('Parent',hpf.fig,'Title','Unit cell parameters','units','normalized','Position',[h_sec v_sec2+br 1-h_sec-br 1-v_sec2-(3*br)],'ShadowColor',[0 0 0],'ForegroundColor',[0 0 0],'HighlightColor',[0.95 0.95 0.95],'BackgroundColor',[0.8 0.8 0.8]);
hpf.inf.tblA = uitable('Parent',hpf.inf.panA,'units','normalized','Position',[br br 1-br 1-br*2]);

% Projected Coordinates Panel
hpf.inf.panP = uipanel('Parent',hpf.fig,'Title','Projected Coordinates','units','normalized','Position',[h_sec v_sec+br 1-h_sec-br v_sec2-v_sec],'ShadowColor',[0 0 0],'ForegroundColor',[0 0 0],'HighlightColor',[0.95 0.95 0.95],'BackgroundColor',[0.8 0.8 0.8]);
hpf.inf.tblP = uitable('Parent',hpf.inf.panP,'units','normalized','Position',[br br 1-br 1-br*4]);


% Parameter Panel & Controls
hpf.par.pan = uipanel('Parent',hpf.fig,'Title','Orientation & Controls','units','normalized','Position',[br br*2 1-br*2 v_sec-br],'ShadowColor',[0 0 0],'ForegroundColor',[0 0 0],'HighlightColor',[0.95 0.95 0.95],'BackgroundColor',[0.8 0.8 0.8]);

btn_width = 1/(5+br*7);
hpf.par.load_cif = uicontrol('Parent',hpf.par.pan,'Style','pushbutton','String','Load Cif file','units','normalized','Position',[br br*6 btn_width 0.3],'FontSize',10);
hpf.par.close = uicontrol('Parent',hpf.par.pan,'Style','pushbutton','String','Cancel','units','normalized','Position',[br+btn_width br*6 btn_width 0.3],'FontSize',10);
hpf.par.save_mat = uicontrol('Parent',hpf.par.pan,'Style','pushbutton','String','Save projected coordinates','units','normalized','Position',[br+btn_width*2 br*6 btn_width 0.3],'FontSize',10);
hpf.par.save_txt = uicontrol('Parent',hpf.par.pan,'Style','pushbutton','String','Export StatSTEM DB file','units','normalized','Position',[br+btn_width*3 br*6 btn_width 0.3],'FontSize',10);
hpf.par.bg = uibuttongroup('Parent',hpf.par.pan,'BorderType','etchedout','Position',[br br+0.4 btn_width 0.3]);
hpf.par.tb1 = uicontrol('Style','radiobutton','String','(uvw)','Parent',hpf.par.bg,'units','normalized','Position',[br br 0.5 1]);
hpf.par.tb2 = uicontrol('Style','radiobutton','String','[hkl]','Parent',hpf.par.bg,'units','normalized','Position',[0.5+br br 0.5 1]);
hpf.par.edit_vec = uicontrol('Style','edit','Parent',hpf.par.pan,'String','0 0 1','units','normalized','Position',[btn_width+br br+0.4 btn_width 0.3]);

% Callbacks
set(hpf.par.close,'Callback',{@closeFig,hpf.fig})
set(hpf.par.save_mat,'Callback',{@save_mat,hpf})
set(hpf.par.save_txt,'Callback',{@save_txt,hpf})
set(hpf.par.load_cif,'Callback',{@load_cif,hpf})
set(hpf.par.edit_vec,'Callback',{@align_vec,hpf})
set(hpf.par.tb1,'Callback',{@align_vec,hpf})
set(hpf.par.tb2,'Callback',{@align_vec,hpf})

function align_vec(~,~,hpf)
    global crystal_par R
    T_hkl = str2num(hpf.par.edit_vec.String); %#ok<ST2NM>
                                         
    [atoms, R, crystal_par] = tfm_align_duplicate_cut(crystal_par, T_hkl, 0, 0, 0, 0, false, hpf.par.tb2.Value);
    
    %Plot projection
    cla(hpf.image.ax2)
    tfm_plot_crystal(atoms, 'g', [R [0 0 0]'], 'h', hpf.image.ax2,'2d');
    fcn_set_proj_coordinates(hpf, atoms)
end

function hpf = load_cif(~, ~, hpf)
    global crystal_par
    [file,cif_path] = uigetfile('*.cif');
    [~,n,~] = fileparts(file);
    if ~isempty(file)
        cla(hpf.image.ax1)
        crystal_par = tfm_get_uc_from_cif([cif_path filesep file]);
        crystal_par.name = n;
        g = tfm_direct_structure_matrix(crystal_par.a, crystal_par.b, crystal_par.c,... 
           crystal_par.alpha, crystal_par.beta, crystal_par.gamma);
       crystal_par.atoms(:,2:4) = crystal_par.atoms(:,2:4)*g; 
       
       tfm_plot_crystal(crystal_par.atoms, 'g', [g [0 0 0]'],'h', hpf.image.ax1)
        
       hpf.inf.tblA.RowName = fieldnames(rmfield(crystal_par,{'atoms','transformations','asym_uc'}));
       hpf.inf.tblA.Data = struct2cell(rmfield(crystal_par,{'atoms','transformations','asym_uc'}));
       align_vec([],[],hpf)
    end
%    hpf.inf.tblA.Position(0) = 0;
end

function fcn_set_proj_coordinates(hpf, atoms)
   hpf.inf.tblP.ColumnName = {'Type', 'X', 'Y', 'Z'};
   hpf.inf.tblP.ColumnWidth = {60,60,60,60};
   hpf.inf.tblP.Data = round(atoms(:,1:4),4);
end

function save_txt(~,~,hpf)
    global crystal_par
    proj_coordinates = hpf.inf.tblP.Data;
    elm = tfm_Z_str(proj_coordinates(:,1));
    xyz = proj_coordinates(:,2:4);
       
    [a, b, c, alpha] = get_projected_lattice();
    [~, id_min] = vector_distances([0 0],xyz);
    xyz = xyz - xyz(id_min,:);

    % xyz2 = tfm_loop_dim([proj_coordinates(:,2) xyz],a,2);
    % xyz2 = tfm_loop_dim(xyz2,b,2);
    % figure(3); clf; 
    % scatter3(xyz2(:,2),xyz2(:,3), xyz2(:,4)); hold on;
    % quiver3(0 ,0 ,0 , a(1) , a(2), a(3)); 
    % quiver3(0 ,0 ,0 , b(1) , b(2), b(3)); 
    % quiver3(0 ,0 ,0 , c(1) , c(2), c(3)); 
    % hold off; xlabel('x'); ylabel('y');
    % axis equal;
    % view([0 0 1]);
    
    xyz_frac = round(xyz/[a' b' c'],4);
    formatSpecXY = '%s	%0.3f	%0.3f';
    formatSpecZ = '%s	%0.3f';
    
    txt = {
        ['Database file for StatSTEM'];
        ['length_a	' num2str(vector_length(a),4)];
        ['length_b	' num2str(vector_length(b),4)];
        ['length_c	' num2str(vector_length(c),4)];
        ['angle_ab	' num2str(alpha,4)];
        '';
        'atoms	x	y';
        };
    for ik = 1:length(elm)
         str = {sprintf(formatSpecXY,elm{ik},xyz_frac(ik,1),xyz_frac(ik,2))};
         txt(end+1)=str;
    end
    txt(end+1)={''};
    txt(end+1)={'zInfo'};
    for ik = 1:length(elm)
         str = {sprintf(formatSpecZ,elm{ik},xyz_frac(ik,3))};
         txt(end+1)=str;
    end
    
    default_name = [replace([crystal_par.name '_' hpf.par.edit_vec.String]," ","") '.txt'];
    [file,path] = uiputfile(default_name);
    if path ~= 0
        writecell(txt,[path filesep file],'Delimiter','tab','QuoteStrings',0)
    end
end

function save_mat(~,~,hpf)
    global crystal_par
    proj_coordinates = hpf.inf.tblP.Data;
    default_name = replace([crystal_par.name '_' hpf.par.edit_vec.String]," ","");
    uisave({'proj_coordinates'},default_name)
end

function closeFig(~,~,fig)
    close(fig)
end    

function [a, b, c, alpha] = get_projected_lattice()
    global R
    % decompose oriented lattice matrix into vectors
    ap = R(1,:);
    bp = R(2,:);
    cp = R(3,:);
    % compute projected areas
    A = [fcn_vec_area_parallelogram(ap,bp),...
        fcn_vec_area_parallelogram(ap,cp),...
    	fcn_vec_area_parallelogram(bp,cp)];

    [~, iA] = max(A);

    % use the vector pair spanning the smallest projected area
    % as new 2D lattice vectors 
    if iA == 1
        a = ap;
        b = bp;
        c = cp;
    elseif  iA == 2
        a = ap;
        b = cp;
        c = bp;
    elseif  iA == 3
        a = bp;
        b = cp;
        c = ap;
    end
    a = [a(1:2) 0];
    b = [b(1:2) 0];
    c = [0 0 sum(c)];
    alpha  = atan2d(norm(cross(a,b)), dot(a,b));
end

function A = fcn_vec_area_parallelogram(a,b)
    ang_ab = fcn_vec_angle_xy(a, b);
    A = vector_length(a(1:2))*vector_length(b(1:2))*sind(ang_ab);
end

function alpha = fcn_vec_angle_xy(a, b)
    a(3) = 0;
    b(3) = 0;
    alpha  = atan2d(norm(cross(a,b)), dot(a,b));
end

function [d, id_min, id_max] = vector_distances(v_ref,xy)
    d = zeros(1,length(xy));
    for id = 1: length(d)
       d(id) = sqrt((v_ref(1)-xy(id,1))^2+(v_ref(2)-xy(id,2))^2);
    end
    [~, id_min] = min(d); 
    [~, id_max] = max(d); 
end

function d = vector_length(v)
    d = sqrt(sum(v.^2));
end