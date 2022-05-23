% Function for Octave compatibility
function b = fcn_startsWith(str, pat)
    for ix = 1:numel(str)
        m = strfind(str{ix}, pat);
        if ~isempty(m) 
            if m(1) == 1
                b(ix) = true;
            else
                b(ix) = false;
            end
        else
            b(ix) = false;
        end
    end
end