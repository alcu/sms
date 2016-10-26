  function x = flipdims(x, varargin)
%|function x = flipdims(x, varargin)
%|
%| Generalization of flipdim() that flips all dimensions.
%| This is useful for the adjoint of convolution operator.
%|
%| option
%|	'odd'	0|1	if 1, make each dimension odd sized after flipping
%|			by appending a zero.  default: 0
%|
%| Jeff Fessler

if nargin < 1, help(mfilename), error(mfilename), end

arg.odd = 0;
arg = vararg_pair(arg, varargin);

if nargin < 1, help(mfilename), error(mfilename), end

for ii = 1:ndims(x)
	x = flipdim(x, ii);

	if arg.odd && ~rem(size(x,ii), 2)
		sz = size(x);
		sz(ii) = 1;
		tmp = zeros(sz, class(x));
		x = cat(ii, x, tmp); % append one zero along dim ii
	end
end
