% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function tfm_plot_crystal(atoms, varargin)
    % Plot the crystal structure
    %
    % Syntax:
    %   tfm_plot_crystal(atoms, varargin)
    %
    % Description:
    %   Plot the crystal structure.
    %
    % Input:
    %   atoms (matrix): The atoms of the crystal.
    %       A structure containing the atomic positions and symbols in the form [symbol (Z), x, y, z].
    %   varargin:
    %       'g': direct structure matrix (used for plotting unit cell bounds)
    %       'h': axes handle
    %       'title': title of the plot
    %       '2d': plot structure in 2d projection
           
    idx_g = find(cellfun(@(x) strcmpi(x, 'g') , varargin),1);
    idx_h = find(cellfun(@(x) strcmpi(x, 'h') , varargin),1);
    idx_title = find(cellfun(@(x) strcmpi(x, 'title') , varargin),1);
    b_2d = find(cellfun(@(x) strcmpi(x, '2d') , varargin),1);
    
    if idx_h
        h = varargin{idx_h+1};
        cla(h);
        delete(findall(gcf,'type','annotation'))
    else
        figure();
        h = gca;
    end
    
    xyz = atoms(:,2:4);
    at_types = unique(atoms(:,1));
    [R, C] = tfm_atomic_radius(at_types);
    Z_str = tfm_Z_str(at_types);   
    
    if ~isempty(b_2d) || length(atoms) > 64
        rng_x = max(max(xyz,[],1)-min(xyz,[],1))+2*max(R);
        markerWidth = get_marker_size(R, rng_x, h)/2;
        plt_fcn = @scatter3flt;
    else
        plt_fcn = @scatter3sph;
        markerWidth = R/2;
    end
    
    hold(h, 'on' );
    for ix = 1:numel(at_types)
       b_at = atoms(:,1) == at_types(ix);
       lg(ix) = plt_fcn(h, xyz(b_at,:), markerWidth(ix), C(ix,:)./255);   
    end
    leg = legend(h, lg, Z_str,'FontSize',24);
    leg.AutoUpdate = 'off';
    hold(h, 'off' );
    
    if idx_g
        gg = varargin{idx_g+1};
        g = gg(1:3,1:3);
        if size(gg,2) == 4 
            tr = gg(:,end)';
        else
            tr = [0 0 0];
        end
        fcn_plot_frame(g,h,tr);
        fcn_plot_tripod(tr-1,h);
        fcn_plot_lat_vec(g,h,tr);
    end
    
    axis(h,'equal')
    if b_2d
        grid(h,'on');
        xlabel(h,'x');
        ylabel(h,'y');
        zlabel(h,'z');
    else
        axis(h,'off');    
    end
    
    if idx_title
        str_title = varargin{idx_title+1};
        if contains(str_title,"_{")
            str_title = join(['$' str_title '$']);
        end   
        annotation('textbox', [0, 0.9, 0.2, 0.1], 'string', str_title,'Interpreter','LaTex','FontSize',24,'LineStyle','none')
    end
    
    drawnow()
end

function h_sc = scatter3flt(h, xyz, markerWidth, C)
    h_sc = scatter3(h, xyz(:,1),xyz(:,2),xyz(:,3), markerWidth, C ,'filled','MarkerEdgeColor','k','LineWidth',1);    
    set(h_sc, 'SizeData', markerWidth.^2)
end

function fcn_plot_frame(g, h, tr)
    basis = [[0 1 0];[0 1 1];[0 0 1];[0 0 0]];
    basis2 = [[0 1 1];[0 0 1];[1 0 1];[1 1 1]];
    frame_1 = basis*g+tr;
    frame_2 = circshift(basis,-1,2)*g+tr;
    frame_3 = circshift(basis,-2,2)*g+tr; 
    frame_4 = basis2*g+tr;          
    frame_5 = circshift(basis2,-1,2)*g+tr;
    frame_6 = circshift(basis2,-2,2)*g+tr; 

    hold(h, 'on');
    plot3(h,frame_1(:,1),frame_1(:,2),frame_1(:,3),'k')
    plot3(h,frame_2(:,1),frame_2(:,2),frame_2(:,3),'k')
    plot3(h,frame_3(:,1),frame_3(:,2),frame_3(:,3),'k')
    plot3(h,frame_4(:,1),frame_4(:,2),frame_4(:,3),'k')
    plot3(h,frame_5(:,1),frame_5(:,2),frame_5(:,3),'k')
    plot3(h,frame_6(:,1),frame_6(:,2),frame_6(:,3),'k')
    hold(h, 'off');
end

function [A,B,C] = fcn_plot_lat_vec(g,h,tr)
    Mg = eye(3) * g;
    A = Mg(1,:);
    B = Mg(2,:);
    C = Mg(3,:);
    hold(h, 'on');
    quiver3(h,tr(1), tr(2) ,tr(3) ,A(1), A(2), A(3), 1,'Color','r','LineWidth',3)
    quiver3(h,tr(1), tr(2) ,tr(3) ,B(1), B(2), B(3), 1,'Color','g','LineWidth',3)
    quiver3(h,tr(1), tr(2) ,tr(3) ,C(1) ,C(2) ,C(3), 1,'Color','b','LineWidth',3)
    hold(h, 'off');
end

function fcn_plot_tripod(offset,h)
    hold(h, 'on');
    quiver3(h,offset(1), offset(2) ,offset(3), 1 ,0 ,0 ,3,'Color','r','LineWidth',3)
    quiver3(h,offset(1), offset(2) ,offset(3), 0 ,1 ,0 ,3,'Color','g','LineWidth',3)
    quiver3(h,offset(1), offset(2) ,offset(3), 0 ,0 ,1 ,3,'Color','b','LineWidth',3)
    hold(h, 'off');
end

function markerWidth = get_marker_size(R, rng_xyz, h)
    currentunits = get(h,'Units');
    set(h, 'Units', 'Points');
    axpos = get(h,'Position');
    set(h, 'Units', currentunits);
    % Calculate Marker width in points
    markerWidth = R./rng_xyz.*axpos(3); 
end