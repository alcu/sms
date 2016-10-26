function [x xinit] = rec_nufft(ksp_dat, kinfo, varargin)
% function x = rec_nufft(ksp_dat, kinfo, varargin)
%
% Reconstructs images from single-shot k-space data using CG with NUFFTs.
%
% Inputs:
% ksp_dat         [num_ksp_data_pts, ncoils, nslcs]
% kinfo           kinfo struct from rec_setup1():
%                 * kinfo.k: [num_ksp_data_pts] Complex, should go from -nx/2 to nx/2.
%                 * kinfo.kdens: [num_ksp_data_pts] Density correction.
%                 * kinfo.kmod: [num_ksp_data_pts] For shifting the final image.
%                 * kinfo.t: [num_ksp_data_pts] Time vector for readout.
% 'scaninfo'      scaninfo struct from rec_setup1().
% 'imsize'        [nx, ny]
% 'chol'          0|1 Use non-iterative cholesky decomposition method, mainly used
%                 for debugging. Default = 0.
% 'nufft'         cell arguments for nufft_init() via Gnufft()
%                 (excluding the first argment, omega)
% 'mask'          [nx, ny, nslcs]
% 'dofmapcorr'    0|1 Do fmap correction or not.
% 'fm'            [nx, ny, nslcs] (Hz). Has no effect if dofmapcorr == 0.
% CG params:
% 'xinit'         [nx, ny, nslcs, ncoils]
% 'beta'          Regularization param.
% 'niter'         Num CG iters.
% 'stop_diff_tol' Tolerance for CG stoppage.
%
% Outputs:
% x               [nx, ny, nslcs, ncoils]

%%% Defaults.
arg.scaninfo.revflg = 1; % 0 = spiral-out (forward), 1 = spiral-in (reverse).
arg.scaninfo.te = 30; % (ms)
arg.scaninfo.ngap1 = 0;
arg.imsize = [64 64];
%! Hack needed b/c other default params depend on these default params.
argtmp = vararg_pair(arg, varargin, 'allow_new', 1);
arg.imsize = argtmp.imsize;
%! End hack.
nx = arg.imsize(1);
ny = arg.imsize(2);
nslcs = size(ksp_dat, 3);
ncoils = size(ksp_dat, 2);
arg.chol = 0;
arg.nufft = {[nx ny], [6 6], 2*[nx ny], [nx/2 ny/2], 'minmax:kb'};
arg.mask = false([nx ny nslcs]);
[xx yy] = ndgrid([-nx/2:nx/2 - 1], [-ny/2:ny/2 - 1]); arg.mask(:,:, 1) = sqrt(xx.^2 + yy.^2) <= nx/2; % Simple mask.
for slc = 2:nslcs; arg.mask(:,:, slc) = arg.mask(:,:, 1); end; % Same mask for all slices.
arg.dofmapcorr = 1;
arg.fm = zeros([nx ny nslcs]);
%arg.dt = 2.5e-6;
% CG params.
arg.xinit = zeros([nx ny nslcs ncoils]);
arg.beta = 273;%241.7618716942761;%400;%1200;
arg.niter = 3000;
arg.stop_diff_tol = 0.0008; % 0 for no stoppage based on tolerance.

%%% Replace defaults with user inputs.
arg = vararg_pair(arg, varargin);

%%% Adjust for reverse (spiral-in) spiral.
if arg.scaninfo.revflg == 1
    tts_shift = arg.scaninfo.te*1e-3 - arg.scaninfo.ngap1*1e-6/2;
    tts = -kinfo.t(1:end) + tts_shift;
elseif arg.scaninfo.revflg == 0
    tts_shift = arg.scaninfo.te*1e-3 + arg.scaninfo.ngap1*1e-6/2;
    tts = kinfo.t(1:end) + tts_shift;
elseif arg.scaninfo.revflg == 3
    tts = kinfo.t;
else
    error(['Not sure what you meant by scaninfo.revflg = ' num2str(arg.scaninfo.revflg)])
end

%ti = [0:1:size(ksp_dat, 1) - 1]'*arg.dt; % k-space sample times (for field map correction).
x = zeros(nx, ny, nslcs, ncoils);
xinit = zeros(nx, ny, nslcs, ncoils);
for slc = 1:nslcs
    %prt ['\nDoing slice ' num2str(slc) ' of ' num2str(nslcs) '...\n']
    prt ['********* Doing slice ' num2str(slc) ' of ' num2str(nslcs)]
    
    % Generate NUFFT object for each slice.
    if arg.dofmapcorr
        if isequal(arg.fm(:,:, slc), zeros(nx, ny))
            fmap_cur = [];
        else
            fmap_cur = arg.fm(:,:, slc);
        end
        
        dft_mat = Gmri([real(kinfo.k.')/nx imag(kinfo.k.')/ny], arg.mask(:,:, slc), ...
                       'fov', [nx ny], 'basis', {'rect'}, ...
                       'nufft', arg.nufft, ...
                       'ti', tts, 'zmap', i*2*pi*fmap_cur, 'L', 6);
        
        if isequal(arg.xinit, zeros(size(arg.xinit)))
            dft_mat_dir = feval(dft_mat.arg.new_image_basis, dft_mat, {'dirac'});
        end
    else
        dft_mat = Gnufft(arg.mask(:,:, slc), {[pi/(nx/2)*real(kinfo.k.') pi/(ny/2)*imag(kinfo.k.')], arg.nufft{:}});
    end
    
    if arg.chol
        R = Reg1(arg.mask(:,:, slc), 'type_penal', 'mat', 'type_diff', 'spmat', 'offsets', penalty_offsets([], [nx ny]), 'order', 1, 'beta', arg.beta, 'distance_power', 1);
    else
        R = Reg1(arg.mask(:,:, slc), 'type_penal', 'def', 'type_diff', 'mex', 'offsets', penalty_offsets([], [nx ny]), 'order', 1, 'beta', arg.beta, 'distance_power', 1);
    end
    
    for coil = 1:ncoils
        prt ['Doing coil ' num2str(coil) ' of ' num2str(ncoils)]
        
        if ~isequal(arg.xinit, zeros(size(arg.xinit))) | ~arg.dofmapcorr
            init = arg.xinit(:,:, slc, coil);
            init = init(arg.mask(:,:, slc));
        else
            init = (dft_mat_dir' * (ksp_dat(:, coil, slc).*kinfo.kdens.')) / (nx*ny);
            xinit(:, :, slc, coil) = embed(init, arg.mask(:,:, slc));
        end
        
        if arg.scaninfo.revflg == 1
            ksp_dat(:, coil, slc) = conj(kinfo.kmod).'.*ksp_dat(end:-1:1, coil, slc);
        else
            ksp_dat(:, coil, slc) = kinfo.kmod.'.*ksp_dat(:, coil, slc);
        end
        
        if arg.chol
            tic
            A = sparse(dft_mat);
            prt ['Gen sparse A: ' num2str(toc)]
            tic
            C = R.C;
            prt ['Gen C: ' num2str(toc)]
            tic
            H = (A' * A + C' * C);
            prt ['Gen H: ' num2str(toc)]
            tic
            %xtmp = full(dft_mat' * dft_mat + R.C' * R.C) \ (dft_mat' * ksp_dat(:, coil, slc));
            xtmp = H \ (dft_mat' * ksp_dat(:, coil, slc));
            prt ['Compute xtmp: ' num2str(toc)]
        else
            [xtmp cginfo] = qpwls_pcg1(init, dft_mat, 1, ksp_dat(:, coil, slc), R.C, 'niter', arg.niter, 'stop_diff_tol', arg.stop_diff_tol);
            prt ['Num iters done for CG: ' num2str(size(cginfo, 1))]
        end
        x(:,:, slc, coil) = embed(xtmp, arg.mask(:,:, slc));
        
    end
end
