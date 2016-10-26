function fn_out = escape_fn(fn)
% function fn_out = escape_fn(fn)
%
% Escapes appropriate characters in fn using bash's printf.
%
% IMPORTANT: Assumes the use of bash as the system's shell. If not using
% bash, you may want to set the MATLAB_SHELL environment variable.
%
% Inputs:
% fn: path string to escape.
%
% Outputs:
% fn_out: escaped path string.

if fn(1) == '~' % Don't escape tildes.
    [status cmdout] = system(['x=''' fn(2:end) '''; x=$(printf ''%q'' "$x"); echo $x']);
    cmdout = ['~' cmdout];
else
    [status cmdout] = system(['x=''' fn '''; x=$(printf ''%q'' "$x"); echo $x']);
end
fn_out = cmdout(1:end - 1); % Last character of cmdout will be a carriage return.
