function ob = ftmodsens_create3(k, mod_kspace, c, fmap, ti, varargin)
% function ob = ftmodsens_create3(k, mod_kspace, c, fmap, ti, varargin)
% Constructs system fatrix2 for SENSE-like simultaneous multislice recon.
%
% Note about mask implementation: See fatrix2.m.
%
% Inputs:
% k               Vector of complex k-space frequencies. (rad/s)
% mod_kspace      [length(k), ncoils, nslcs]
% c               [nx, ny, nslcs, ncoils] Coil sensitivities.
% fmap            [nx, ny, nslcs]
% ti              [length(k)] Time vector during readout. (s)
% 'mask'          [nx, ny, nslcs]
% 'nufft'         cell arguments for nufft_init() via Gnufft()
%                 (excluding the first argment, omega)

%sizec = size(c);
nx = size(c, 1);
ny = size(c, 2);
arg.nslcs = size(c, 3);
arg.ncoils = size(c, 4);

%%% Defaults.
arg.mask = true(nx, ny, arg.nslcs);
arg.nufft = {[nx ny], [6 6], 2*[nx ny], [nx/2 ny/2], 'minmax:kb'};

%%% Replace defaults with user inputs.
arg = vararg_pair(arg, varargin);

arg.dft_mat = cell(arg.nslcs, 1);
arg.nksamp = length(k);
%dt = 4e-6; % Depends on bandwidth of scanner data sampling (sampling rate).
%ti = [0:1:arg.nksamp - 1]'*dt; % k-space sample times (for field map correction).
for slc = 1:arg.nslcs
    % arg.dft_mat{slc} = Gnufft(arg.mask(:,:,slc), {[real(k) imag(k)], arg.nufft{:}});
    if isequal(fmap(:,:, slc), zeros(nx, ny))
        fmap_cur = [];
    else
        fmap_cur = fmap(:,:, slc);
    end
    
    arg.dft_mat{slc} = Gmri([real(k) imag(k)]/(2*pi), arg.mask(:,:,slc), ...
                            'fov', size(arg.mask(:,:,slc)), 'basis', {'rect'}, ...
                            'nufft', arg.nufft, ...
                            'ti', ti, 'zmap', i*2*pi*fmap_cur, 'L', 6);
end
arg.mod_kspace = mod_kspace;
arg.c = c;

ob = fatrix2('odim', [arg.nksamp arg.ncoils], 'idim', size(arg.mask), 'mask', arg.mask, 'arg', arg, ...
             'forw', @Atimes, 'back', @Astartimes);

end % function ob = ftmodsens_create3()


function b = Atimes(arg, x)

b = zeros(arg.nksamp, arg.ncoils);
for coil = 1:arg.ncoils
    for slc = 1:arg.nslcs
        b(:, coil) = b(:, coil) + ...
            arg.mod_kspace(:, coil, slc) .* ...
            (arg.dft_mat{slc} * (arg.c(:,:, slc, coil) .* x(:,:, slc)));
    end
end

end % function b = Atimes(arg, x)


function x = Astartimes(arg, b)

x = zeros(size(arg.mask));
for slc = 1:arg.nslcs
    for coil = 1:arg.ncoils
        x(:,:, slc) = x(:,:, slc) + ...
            conj(arg.c(:,:, slc, coil)) .* ...
            embed(arg.dft_mat{slc}' * (conj(arg.mod_kspace(:, coil, slc)) .* b(:, coil)), arg.mask(:,:, slc));
    end
end

end % function x = Astartimes(arg, b)
