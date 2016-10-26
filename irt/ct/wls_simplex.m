 function x = wls_simplex(A, y, Wh, x, varargin)
%function x = wls_simplex(A, y, Wh, x, [options])
%| min_x || Wh * (A x - y) ||^2 + reg || x ||^2
%| subject to simplex constraint: 0 <= x <= 1 and sum(x) = 1
%|
%| based on:
%| x = lsqlin(C,d,A,b,Aeq,beq) solves the least-squares
%| (with equality constraints) problem:
%| min_x || C*x-d ||^2 subject to A*x <= b and Aeq*x = beq
%| x = lsqlin(C,d,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
%| bounds on the design variables, x, so that LB <= x <= UB.
%|
%| option
%|	'inprodv'	if (row vector) provided, require 1 = inprodv * x
%|	'maxiter'	max # iterations for lsqlin (default: 400)
%|	'optimset' {}	cell arguments for optimset().  default: {}
%|	'reg'		regularization parameter (default 1e-6)
%|	'how'		'lsqlin' (default if present), 'pd' (else default)
%|			for primal dual, or 'gp' (gradient projection)
%|			and pd and gp are still under development (todo: test)
%|
%| todo: try keerthi:02:cog
%|
%| Copyright 2006-1-1, Jeff Fessler, University of Michigan

if nargin == 1 && streq(A, 'test'), wls_simplex_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if ~isvar('Wh') || isempty(Wh)
	Wh = 1;
end
if ~isvar('x'), x = []; end

[arg x] = wls_simplex_arg(A, y, Wh, x, varargin{:});

if isempty(arg.how)
	if isfield(optimset, 'LargeScale')
		arg.how = 'lsqlin';
	else
		arg.how = 'pd';
	end
end

switch(arg.how)
case 'lsqlin'
	x = wls_simplex_lsqlin(A, y, Wh, x, varargin{:});
case 'gp'
	x = wls_simplex_gp(A, y, Wh, x, varargin{:});
case 'pd'
	x = wls_simplex_pd(A, y, Wh, x, varargin{:});
otherwise
	fail('method "%s" unknown', arg.how)
end

% check for significant negative values
if any(x < -eps)
	minmax(x)
	warn('zeroing %d of %d "big" negatives', sum(x < -eps), numel(x))
	x = max(x,0);
	x = x / sum(x);
end

% purge negligible negative values
if any(x < 0)
	printm('zeroing %d of %d tiny negatives', sum(x < 0), numel(x))
	x = max(x,0);
	x = x / sum(x);
end

if any(x > 1)
	fail 'bug: x > 1'
end

x = x / sum(x);


% wls_simplex_arg()
function [arg x] = wls_simplex_arg(A, y, Wh, x, varargin)

arg.how = [];
arg.inprodv = [];
arg.maxiter = 400;
arg.optimset = {};
arg.reg = 1e-16;
arg.tol = 1e-7;
arg = vararg_pair(arg, varargin);

[m n] = size(A);

if ~isequal([m 1], size(y))
	fail 'y must be [nrow(A) 1]'
end

if ~isscalar(Wh) && ~isequal([m m], size(Wh))
	fail 'Wh must be [nrow(A) nrow(A)]'
end

if isempty(x)
	x = ones(n,1) / n;
else
	if ~isequal([n 1], size(x))
		fail 'x must be [n 1]'
	end
	x = max(x, 0);
	x = x / sum(x);
end


% wls_simplex_pd()
% primal-dual method based on chambolle and pock
function x = wls_simplex_pd(A, y, Wh, x, varargin)
arg = wls_simplex_arg(A, y, Wh, x, varargin{:});
tol = arg.tol;

B = Wh * A;
bb = B' * (Wh * y);
hess = B' * B + arg.reg * eye(ncol(A));
maxeig = max(eig(hess));
if ~isempty(arg.inprodv)
	fail 'inprodv not done'
end

step_factor = 0.99;
sigma = step_factor / sqrt(maxeig);
tau = sigma;
theta = 1.0; % todo

if isscalar(Wh)
	mat1 = Wh^2 / (sigma + Wh^2); % inv(sig I + W) * W
%elseif isequal(diag(diag(Wh)), Wh)
else
	keyboard
	fail 'todo'
end

%A = [Wh * A; sqrt(arg.reg) * eye(ncol(A))];
%y = [Wh * y; zeros(ncol(A),1)];

z = 0;
xb = x;

for iter=1:arg.maxiter
	% | Wh * (A x - y) |^2 + reg | x |^2
	old = x;
	z = mat1 * (z + sigma * (A * xb - y));
	x = x - tau * (A' * z);
	x = ir_project_simplex(x / (1 + tau * arg.reg)); 
	xb = x + theta * (x - old);
	if 0 && ~mod(iter, 1000)
		plot(x), axisy([0 1]), drawnow
	end
	if norm(x - old) / norm(x) < tol
		printm('stopped at %d iterations for tol=%g', iter, tol)
		return
	end
end
printm('reached maximum %d iterations with change=%g', ...
	iter, norm(x - old) / norm(x))


% wls_simplex_gp()
% gradient projection method
function x = wls_simplex_gp(A, y, Wh, x, varargin)
arg = wls_simplex_arg(A, y, Wh, x, varargin{:});
tol = arg.tol;

B = Wh * A;
bb = B' * (Wh * y);
hess = B' * B + arg.reg * eye(ncol(A));
step = 1 / max(eig(hess));
if ~isempty(arg.inprodv)
	fail 'inprodv not done'
end

for iter=1:arg.maxiter
	% | Wh * (A x - y) |^2 + reg | x |^2
%	grad = A' * Wh' * Wh * (A * x - y) + reg * x;
	grad = hess * x - bb;
	old = x;
	x = x - step * grad;
	x = ir_project_simplex(x); 
	if 0 && ~mod(iter,10)
		plot(x), axisy([0 1]), drawnow
	end
	if norm(x - old) / norm(x) < tol
		printm('stopped at %d iterations for tol=%g', iter, tol)
		return
	end
end
printm('reached maximum %d iterations with relerr=%g', ...
	iter, norm(x - old) / norm(x))


% wls_simplex_lsqlin()
function x = wls_simplex_lsqlin(A, y, Wh, x, varargin)
arg = wls_simplex_arg(A, y, Wh, x, varargin{:});

n = ncol(A);

if isempty(arg.inprodv)
	Aeq = ones(1,n);
	Beq = 1;
else
	if ~isequal([1 n], size(arg.inprodv))
		fail 'inprodv must be [1 n]'
	end

	if all(arg.inprodv > 1)
		minmax(arg.inprodv)
		fail '<inprodv, x> = 1 infeasibile'
	end

	Aeq = [ones(1,n); arg.inprodv];
	Beq = [1; 1];
end

opt = optimset('largescale', 'off', 'display', 'off', ...
	'maxiter', arg.maxiter, arg.optimset{:});

C = Wh * A;
d = Wh * y;
if arg.reg % trick: build regularizer into matrix
	C = [C; sqrt(arg.reg) * eye(n)];
	d = [d; zeros(n,1)];
end
[x resnorm residual exitflag output lambda] = ...
	lsqlin(C, d, [], [], Aeq, Beq, zeros(1,n), ones(1,n), x, opt);

if exitflag == 0
	warn('lsqlin exitflag=%d, may need more iterations than %d', ...
		exitflag, arg.maxiter)
	if 1
		pr resnorm
		minmax(residual)
		pr output
		pr output.message
		pr lambda
		minmax(lambda.lower)
		minmax(lambda.upper)
	end
elseif exitflag ~= 1
	printm('lsqlin exitflag=%d', exitflag)
	keyboard
end

if isempty(resnorm)
	fail 'inconsistent input bounds'
end


% wls_simplex_test
function wls_simplex_test
ee = linspace(1, 6, 31);
tt = linspace(0, 5, 101)';
A = exp(-tt * ee); % [nt ne]
% y0 = exp(-tt) .* (1 - exp(-3*tt));
x0 = zeros(numel(ee), 1);
x0([7 13 29]) = [0.1 0.6 0.3];
%x0(13) = 1;
y0 = A * x0;
reg = 1e-9;

hows = {'gp', 'pd'};
for ih=2:numel(hows)
	how = hows{ih};
	x1 = wls_simplex(A, y0, [], [], 'how', how, 'reg', reg, 'maxiter', 1e5);
	y1 = A * x1;

	if im
		im plc 2 1
		im subplot 1
		plot(tt, y0, '-', tt, y1, '--')
		im subplot 2
		plot([x0 x1], '-o')
		max_percent_diff(y0, y1)
	end
end

if isfield(optimset, 'LargeScale')
	x2 = wls_simplex(A, y0, [], [], 'how', 'lsqlin', 'reg', reg);
	y2 = A * x2;
	max_percent_diff(y0, y2)
	equivs(x2, x1, 'thresh', 6e-6, 'fail', 0)
	equivs(y2, y1, 'thresh', 3e-6, 'fail', 0)
	if im
		im plc 2 2
		im subplot 1
		plot(tt, [y0 y1 y2]), legend('y0', 'y1', 'y2')
		im subplot 2
		plot(tt, [y1-y0 y2-y0]), legend('y1-y0', 'y2-y0')
		im subplot 3
		plot([x0 x1 x2], '-o'), legend('x0', 'x1', 'x2')
		im subplot 4
		plot([x1-x0 x2-x0]), legend('x1-x0', 'x2-x0')
	end
end
