% Copyright (C) 2022 Thomas Friedrich
% University of Antwerp - All Rights Reserved. 
% You may use, distribute and modify
% this code under the terms of the GPL3 license.
% You should have received a copy of the GPL3 license with
% this file. If not, please visit: 
% https://www.gnu.org/licenses/gpl-3.0.en.html

function [atoms, a, b, c, alpha, beta, gamma, sg, hmg, transformations, formula] = tfm_import_cif(filename)

    str = fileread(filename);
    str = strtrim(splitlines(str));
    
    [~, a] = fcn_find_prm(str,'_cell_length_a');
    [~, b] = fcn_find_prm(str,'_cell_length_b');
    [~, c] = fcn_find_prm(str,'_cell_length_c');
    [~, alpha] = fcn_find_prm(str,'_cell_angle_alpha');
    [~, beta] = fcn_find_prm(str,'_cell_angle_beta');
    [~, gamma] = fcn_find_prm(str,'_cell_angle_gamma');
    [~, sg] = fcn_find_prm(str,'_space_group_it_number');
    if isempty(sg)
      [~, sg] = fcn_find_prm(str,'_symmetry_int_tables_number');
    end 
    [~, formula] = fcn_find_prm(str,'_chemical_name_systematic',true);
    if isempty(formula)
        [~, formula] = fcn_find_prm(str,'_chemical_formula_sum',true);
    end
    [~, hmg] = fcn_find_prm(str,'_symmetry_space_group_name_h-m', true);
    transformations = fcn_read_loops(str,{'_symmetry_equiv_pos_as_xyz'},@fcn_parse_sym_ops);
    atoms = fcn_read_loops(str, {'_atom_site_type_symbol','_atom_site_fract_x','_atom_site_fract_y','_atom_site_fract_z'},@fcn_parse_atoms);
end

function [id, prm] = fcn_find_prm(str, pat, b_str)
    id = fcn_contains(str,pat);
    if sum(id) > 0
        str_cell = strsplit(str{id},{' ','\t'});
        if nargin > 2 && b_str
            prm = replace(strjoin(str_cell(2:end)),"'","");
        else
            prm = fcn_erase_brackets(str_cell{end});
        end   
    else
        prm = [];
    end
    id = find(id);
end

function dat = fcn_read_loops(str, loop_label, parse_fcn)
    id = find(fcn_contains(str,'loop_'));
    n_id = numel(id);
    col = zeros(1,numel(loop_label));
    for il = 1:n_id
        loop_start = id(il);
        tline = strtrim(str{loop_start+1});
        iv = 1;
        while fcn_startsWith(tline,'_')
           va{iv} = tline;
           iv = iv + 1;
           tline = strtrim(str{loop_start + iv}); 
        end
        b_at = fcn_find_prm(va,loop_label{1});
        if ~isempty(b_at)
           for c = 1:numel(loop_label)
                col(c) = fcn_find_prm(va,loop_label{c});
           end
           ir = 1;
           while true && loop_start + iv +ir-1 <= numel(str)
               tline = strtrim(str{loop_start + iv +ir-1});
               if fcn_startsWith(tline,'_') || fcn_startsWith(tline,'loop_') || isempty(strtrim(tline))
                   break
               end
               tline = replace(tline,'"',"'");
               sub_str = regexp(tline,"'(.*?)'",'match');
               if ~isempty(sub_str)
                   for s = 1: numel(sub_str)
                        tline = replace(tline,sub_str{s},replace(sub_str{s},{' ',''''},''));
                   end
               end
               cell_arr = strsplit(tline,{' ','\t'});
               dat(ir,:) = parse_fcn(cell_arr(col));
               ir = ir + 1;
           end
        end
    end
end

function atom = fcn_parse_atoms(cell_arr)
    typ = char(cell_arr(1));
    atom(1) = ilm_Z(typ(isletter(typ)));
    atom(2:4) = cellfun(@fcn_erase_brackets, cell_arr(2:4), 'UniformOutput', true);
end

function sym_op = fcn_parse_sym_ops(cell_arr)
    sym_op = split(string(cell_arr),',');
end

% Octave compatibility functions
function flt = fcn_erase_brackets(str)
    [s,e] = regexp(str,'(\(\d+\))');
    if ~isempty(s)
      str(s:e) = [];
    end
    flt = str2double(str);
end

function b = fcn_contains(str,pat)
  for ix = 1:numel(str)
    m = strfind(lower(str{ix}), lower(pat));
    b(ix) = ~isempty(m); 
  end
end

function b = fcn_startsWith(str,pat)
  m = strfind(str, pat);
  if ~isempty(m) 
    if m(1) == 1
      b = true;
    else
      b = false;
    end
  else
    b = false;
  end
end

function s = splitlines(string)
  s = strsplit(string, '\n');
  s = reshape(s,size(s, 2), size(s,1));
end