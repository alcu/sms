function fm = rec_fmap(ksp_dat0, ksp_dat1, kinfo, varargin)
% function fm = rec_fmap(ksp_dat0, ksp_dat1, kinfo, varargin)
%
% Constructs field maps from single-shot k-space data using CG with NUFFTs.
%
% Inputs:
% ksp_dat0        [num_ksp_data_pts, ncoils, nslcs]
% ksp_dat1        [num_ksp_data_pts, ncoils, nslcs] K-space data for second volume, acquired with a different TE.
% kinfo           kinfo struct from rec_setup1():
%                 * kinfo.k: [num_ksp_data_pts] Complex, should go from -nx/2 to nx/2.
%                 * kinfo.kdens: [num_ksp_data_pts] Density correction (not necessary).
%                 * kinfo.kmod: [num_ksp_data_pts] For shifting the final image.
%                 * kinfo.t: [num_ksp_data_pts] Time vector for readout.
% 'scaninfo'      scaninfo struct from rec_setup1().
% 'args'          args struct from rec_setup1().
% 'imsize'        [nx, ny]
% 'nufft'         cell arguments for nufft_init() via Gnufft()
%                 (excluding the first argment, omega)
% 'mask'          [nx, ny, nslcs]
% 'fm_reps'       Number of field map construction repetitions.
% 'fm'            [nx, ny, nslcs] (Hz). Initial fmap guess.
% CG params:
% 'xinit'         [nx, ny, nslcs, ncoils]
% 'beta'          Regularization param.
% 'niter'         Num CG iters.
% 'stop_diff_tol' Tolerance for CG stoppage.
%
% Outputs:
% fm              [nx, ny, nslcs]

%%% Defaults.
arg.scaninfo.revflg = 1; % 0 = spiral-out (forward), 1 = spiral-in (reverse).
arg.scaninfo.te = 30; % (ms)
arg.scaninfo.ngap1 = 0;
arg.scaninfo.mapdel = 2000;
arg.args.fsize = 7;
arg.imsize = [64 64];
%! Hack needed b/c other default params depend on these default params.
argtmp = vararg_pair(arg, varargin, 'allow_new', 1);
arg.imsize = argtmp.imsize;
%! End hack.
nx = arg.imsize(1);
ny = arg.imsize(2);
nslcs = size(ksp_dat0, 3);
ncoils = size(ksp_dat0, 2);
arg.nufft = {[nx ny], [6 6], 2*[nx ny], [nx/2 ny/2], 'minmax:kb'};
arg.mask = false([nx ny nslcs]);
[xx yy] = ndgrid([-nx/2:nx/2 - 1], [-ny/2:ny/2 - 1]); arg.mask(:,:, 1) = sqrt(xx.^2 + yy.^2) <= nx/2; % Simple mask.
for slc = 2:nslcs; arg.mask(:,:, slc) = arg.mask(:,:, 1); end; % Same mask for all slices.
arg.fm_reps = 2;
arg.fm = zeros([nx ny nslcs]);
%arg.dt = 2.5e-6;
% CG params.
arg.xinit = zeros([nx ny nslcs ncoils]);
arg.beta = 273;%241.7618716942761;%400;%1200;
arg.niter = 3000;
arg.stop_diff_tol = 0.0008; % 0 for no stoppage based on tolerance.

%%% Replace defaults with user inputs.
arg = vararg_pair(arg, varargin);

for ii = 1:arg.fm_reps
    prt char(['********* Doing fmap rep ' num2str(ii) ' of ' num2str(arg.fm_reps) ' *********'])
    
    prt char('********* Recon vol 0: *********')
    imarall0 = rec_nufft(ksp_dat0, kinfo, 'scaninfo', arg.scaninfo, 'imsize', arg.imsize, 'nufft', arg.nufft, 'mask', arg.mask, 'fm', arg.fm, 'xinit', arg.xinit, 'beta', arg.beta, 'niter', arg.niter, 'stop_diff_tol', arg.stop_diff_tol);
    prt char('********* Recon vol 1: *********')
    imarall1 = rec_nufft(ksp_dat1, kinfo, 'scaninfo', arg.scaninfo, 'imsize', arg.imsize, 'nufft', arg.nufft, 'mask', arg.mask, 'fm', arg.fm, 'xinit', arg.xinit, 'beta', arg.beta, 'niter', arg.niter, 'stop_diff_tol', arg.stop_diff_tol);

    arg.fm = sum(imarall1 .* conj(imarall0), 4); % Size: [nx ny nslcs]

    % Add fancy smoother.
    mnfm = mean(reshape(arg.fm, nx*ny, nslcs), 1); % Size: [1 nslcs]
    fmtmp = reshape(angle(bsxfun(@times, reshape(arg.fm, nx*ny, nslcs), exp(-i*angle(mnfm)))), nx, ny, nslcs); % Size: [nx ny nslcs]
    fmmag = sqrt(abs(arg.fm)); % Size: [nx ny nslcs]
    fmsm = zeros(nx, ny, nslcs);
    for slc = 1:nslcs
        fmsm(:,:, slc) = angle(conv2(fmmag(:,:, slc).*exp(i*fmtmp(:,:, slc)), ones(arg.args.fsize), 'same')) + angle(mnfm(slc)); % Size: [nx ny nslcs]
    end
    arg.fm = fmsm/(arg.scaninfo.mapdel*1e-6)/2/pi; % in Hz

end
fm = arg.fm;
