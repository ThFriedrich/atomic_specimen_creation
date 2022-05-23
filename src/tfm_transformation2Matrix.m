% https://github.com/gudasergey/fdmnes/blob/master/source/spgroup.f90
% Based on the fdmnes code by Sergey Guda (no license specified at this
% time)

function [R,t] = tfm_transformation2Matrix(tr)
    R = zeros(3);
    t = zeros(1,3);
    for i = 1:3
        switch true 
            case any(strcmpi(tr(i),{'x','+x'}))
                R(i,1) = 1;

            case any(strcmpi(tr(i),{'-x'}))
                R(i,1) = -1;

            case any(strcmpi(tr(i),{'y','+y'}))
                R(i,2) = 1;

            case any(strcmpi(tr(i),{'-y'}))
                R(i,2) = -1;

            case any(strcmpi(tr(i),{'z','+z'}))
                R(i,3) = 1;

            case any(strcmpi(tr(i),{'-z'}))
                R(i,3) = -1;

            case any(strcmpi(tr(i),{'x-y','+x-y','-y+x'}))
                R(i,1) = 1;
                R(i,2) = -1;

            case any(strcmpi(tr(i),{'y-x','+y-x','-x+y'}))
                R(i,1) = -1;
                R(i,2) = 1;

            case any(strcmpi(tr(i),{'1/2+x','x+1/2','+x+1/2'}))
                R(i,1) = 1;
                t(i) = 0.5;

            case any(strcmpi(tr(i),{'x-1/2','+x-1/2','-1/2+x'}))
                R(i,1) = 1;
                t(i) = -0.5;

            case any(strcmpi(tr(i),{'1/2-x','-x+1/2'}))
                R(i,1) = -1;
                t(i) = 0.5;

            case any(strcmpi(tr(i),{'-1/2-x','-x-1/2'}))
                R(i,1) = -1;
                t(i) = -0.5;

            case any(strcmpi(tr(i),{'1/2+y','y+1/2','+y+1/2'}))
                R(i,2) = 1;
                t(i) = 0.5;

            case any(strcmpi(tr(i),{'1/2-y','-y+1/2'}))
                R(i,2) = -1;
                t(i) = 0.5;

            case any(strcmpi(tr(i),{'-1/2+y','y-1/2','+y-1/2'}))
                R(i,2) =  1;
                t(i) = -0.5;

            case any(strcmpi(tr(i),{'-y-1/2','-1/2-y'}))
                R(i,2) = -1;
                t(i) = -0.5;

            case any(strcmpi(tr(i),{'z+1/2','+z+1/2','1/2+z'}))
                R(i,3) = 1;
                t(i) = 0.5;

            case any(strcmpi(tr(i),{'1/2-z','-z+1/2'}))
                R(i,3) = -1;
                t(i) = 0.5;

            case any(strcmpi(tr(i),{'-1/2+z','z-1/2','+z-1/2'}))
                R(i,3) = 1;
                t(i) = -0.5;

            case any(strcmpi(tr(i),{'-z-1/2','-1/2-z'}))
                R(i,3) = -1;
                t(i) = -0.5;

            case any(strcmpi(tr(i),{'1/4+x','x+1/4','+x+1/4'}))
                R(i,1) = 1;
                t(i) = 0.25;

            case any(strcmpi(tr(i),{'1/4-x','-x+1/4'}))
                R(i,1) = -1;
                t(i) = 0.25;

            case any(strcmpi(tr(i),{'1/4+y','y+1/4','+y+1/4'}))
                R(i,2) = 1;
                t(i) = 0.25;

            case any(strcmpi(tr(i),{'1/4-y','-y+1/4'}))
                R(i,2) = -1;
                t(i) = 0.25;

            case any(strcmpi(tr(i),{'1/4+z','z+1/4','+z+1/4'}))
                R(i,3) = 1;
                t(i) = 0.25;

            case any(strcmpi(tr(i),{'1/4-z','+1/4-z','-z+1/4'}))
                R(i,3) = -1;
                t(i) = 0.25;

            case any(strcmpi(tr(i),{'3/4+x','x+3/4','+x+3/4'}))
                R(i,1) = 1;
                t(i) = 0.75;

            case any(strcmpi(tr(i),{'3/4-x','+3/4-x','-x+3/4'}))
                R(i,1) = -1;
                t(i) = 0.75;

            case any(strcmpi(tr(i),{'3/4+y','+3/4+y','y+3/4','+y+3/4'}))
                R(i,2) = 1;
                t(i) = 0.75;

            case any(strcmpi(tr(i),{'3/4-y','+3/4-y','-y+3/4'}))
                R(i,2) = -1;
                t(i) = 0.75;

            case any(strcmpi(tr(i),{'3/4+z','+3/4+z','z+3/4','+z+3/4'}))
                R(i,3) = 1;
                t(i) = 0.75;

            case any(strcmpi(tr(i),{'3/4-z','+3/4-z','-z+3/4'}))
                R(i,3) = -1;
                t(i) = 0.75;

            case any(strcmpi(tr(i),{'z+1/6','+z+1/6','1/6+z'}))
                R(i,3) = 1;
                t(i) = 1/6;

            case any(strcmpi(tr(i),{'-z+1/6','+1/6-z','1/6-z'}))
                R(i,3) = -1;
                t(i) = 1/6;

            case any(strcmpi(tr(i),{'z+1/3','+z+1/3','1/3+z'}))
                R(i,3) = 1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'-z+1/3','1/3-z','+1/3-z'}))
                R(i,3) = -1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'-z+5/6','5/6-z','+5/6-z'}))
                R(i,3) = -1;
                t(i) = 5/6;

            case any(strcmpi(tr(i),{'z+2/3','+z+2/3','2/3+z'}))
                R(i,3) = 1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'-z+2/3','2/3-z'}))
                R(i,3) = -1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'z+5/6','+z+5/6','5/6+z','+5/6+z'}))
                R(i,3) = 1;
                t(i) = 5/6;

            case any(strcmpi(tr(i),{'x+1/3','+x+1/3','1/3+x','+1/3+x'}))
                R(i,1) = 1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'x+2/3','+x+2/3','2/3+x'}))
                R(i,1) = 1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'-x+1/3','1/3-x'}))
                R(i,1) = -1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'-x+2/3','2/3-x'}))
                R(i,1) = -1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'y+1/3','+y+1/3','1/3+y'}))
                R(i,2) = 1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'y+2/3','+y+2/3','2/3+y'}))
                R(i,2) = 1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'-y+1/3','1/3-y'}))
                R(i,2) = -1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'-y+2/3','2/3-y'}))
                R(i,2) = -1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'x-y+1/3','+x-y+1/3','-y+x+1/3','1/3+x-y'}))
                R(i,1) = 1;
                R(i,2) = -1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'x-y+2/3','+x-y+2/3','-y+x+2/3','2/3+x-y'}))
                R(i,1) = 1;
                R(i,2) = -1;
                t(i) = 2/3;

            case any(strcmpi(tr(i),{'-x+y+1/3','y-x+1/3','+y-x+1/3','1/3-x+y'}))
                R(i,1) = -1;
                R(i,2) = 1;
                t(i) = 1/3;

            case any(strcmpi(tr(i),{'-x+y+2/3','y-x+2/3','+y-x+2/3','2/3-x+y'}))
                R(i,1) = -1;
                R(i,2) = 1;
                t(i) = 2/3;
        end
    end
end