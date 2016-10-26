function out = flipvec(x)
% function out = flipvec(x)
%
% Simple flip vector function. Uses the appropriate flipdim() call based on
% whether x is a column or row vector.
%
% Obviously, only for vectors.

if isvector(x)
    if size(x, 1) == 1
        out = flipdim(x, 2);
    else
        out = flipdim(x, 1);
    end
else
    error('flipvec() only works on vectors.')
end
