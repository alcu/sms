function imarev = rec_sense(ksp_dat, kinfo, c, varargin)
% function imarev = rec_sense(ksp_dat, kinfo, c, varargin)
%
% Reconstructs images from single-shot k-space data using SENSE (CG with NUFFTs).
% Can also reconstruct simultaneous multislice data.
%
% Inputs:
% ksp_dat         [num_ksp_data_pts, ncoils, nacqs]
% kinfo           kinfo struct from rec_setup1():
%                 * kinfo.k: [num_ksp_data_pts] Complex, should go from -nx/2 to nx/2.
%                 * kinfo.kdens: [num_ksp_data_pts] Density correction (not necessary).
%                 * kinfo.kmod: [num_ksp_data_pts] For shifting the final image.
%                 * kinfo.t: [num_ksp_data_pts] Time vector for readout.
% c               [nx, ny, nslcs_overall, ncoils, nmaps] Sensitivity maps. (nmaps
%                 is the number of sets of ncoils maps produced by ESPIRiT.)
% 'scaninfo'      scaninfo struct from rec_setup1().
% 'imsize'        [nx, ny]
% 'nufft'         cell arguments for nufft_init() via Gnufft()
%                 (excluding the first argment, omega)
% 'mask'          [nx, ny, nslcs_overall, nmaps]
% 'fm'            [nx, ny, nslcs_overall] (Hz).
% 'gzmod'         [num_ksp_data_pts, ncoils, nslcs_overall/nslcs == nacqs, nslcs] (rad)
% CG params:
% 'xinit'         [nx, ny, nslcs_overall]
% 'beta'          Regularization param.
% 'niter'         Num CG iters.
% 'stop_diff_tol' Tolerance for CG stoppage.
%
% Outputs:
% imarev          [nx, ny, nslcs_overall, nmaps]

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
ndat = size(ksp_dat, 1);
ncoils = size(ksp_dat, 2);
nacqs = size(ksp_dat, 3); % Number of data acquisitions for this volume.
nslcs_overall = size(c, 3); % Overall number of slices (after separation).
nslcs = nslcs_overall/nacqs; % Number of simultaneously acquired slices.
arg.nmaps = size(c, 5);
arg.nufft = {[nx ny], [6 6], 2*[nx ny], [nx/2 ny/2], 'minmax:kb'};
arg.mask = false([nx ny nslcs_overall]);
[xx yy] = ndgrid([-nx/2:nx/2 - 1], [-ny/2:ny/2 - 1]); arg.mask(:,:, 1) = sqrt(xx.^2 + yy.^2) <= nx/2; % Simple mask.
for slc = 2:nslcs_overall; arg.mask(:,:, slc) = arg.mask(:,:, 1); end; % Same mask for all slices.
for map = 1:arg.nmaps; arg.mask(:,:,:, map) = arg.mask(:,:,:, 1); end; % Same mask for all components.
arg.fm = zeros([nx ny nslcs_overall]);
%arg.dt = 2.5e-6;
arg.gzmod = zeros(ndat, ncoils, nacqs, nslcs); % (rad)
% CG params.
arg.xinit = zeros([nx ny nslcs_overall]);
arg.beta = 273;%241.7618716942761;%400;%1200;
arg.niter = 3000;
arg.stop_diff_tol = 0.0008; % 0 for no stoppage based on tolerance.
plot_extra = 0;

%%% Replace defaults with user inputs.
arg = vararg_pair(arg, varargin);

%%% Adjust for reverse (spiral-in) spiral.
if arg.scaninfo.revflg == 1
    tts_shift = arg.scaninfo.te*1e-3 - arg.scaninfo.ngap1*1e-6/2;
    tts = -kinfo.t(1:end) + tts_shift;
    arg.gzmod = arg.gzmod(end:-1:1, :, :, :);
    for acq = 1:nacqs
        for coil = 1:ncoils
            ksp_dat(:, coil, acq) = conj(kinfo.kmod).'.*ksp_dat(end:-1:1, coil, acq);
        end
    end
elseif arg.scaninfo.revflg == 0
    tts_shift = arg.scaninfo.te*1e-3 + arg.scaninfo.ngap1*1e-6/2;
    tts = kinfo.t(1:end) + tts_shift;
    for acq = 1:nacqs
        for coil = 1:ncoils
            ksp_dat(:, coil, acq) = kinfo.kmod.'.*ksp_dat(:, coil, acq);
        end
    end
elseif arg.scaninfo.revflg == 3
    tts = kinfo.t;
    for acq = 1:nacqs
        for coil = 1:ncoils
            ksp_dat(:, coil, acq) = kinfo.kmod.'.*ksp_dat(:, coil, acq);
        end
    end
else
    error(['Not sure what you meant by scaninfo.revflg = ' num2str(arg.scaninfo.revflg)])
end

imarev = zeros(nx, ny, nslcs_overall, arg.nmaps);
dolcurves = 0;
for acq = 1:nacqs
    %{
    dolcurves = 1;
    betas = 10.^[-3 : 0.3 : 3.0];
    lcur = zeros(length(betas), 2); % Format: lcur(ii, [moderr prerr])
    for ii = 1:length(betas)
    arg.beta = betas(ii);
    prt ['********* Doing beta number ' num2str(ii) ' of ' num2str(length(betas))]
    %}
    
    prt ['********* Doing acquisition ' num2str(acq) ' of ' num2str(nacqs)]
    
    % Actual slices we are reconstructing.
    slcs_actual = [1 : nacqs : 1 + nacqs*(nslcs - 1)] + acq - 1;
    
    % Use the appropriate experimentally acquired modulation waveform.
    mod_kspace = exp(i*squeeze(arg.gzmod(:, :, acq, :)));
    
    %%% Create fatrix AC (system matrix for Fessler's CG function: qpwls_pcg).
    AC = ftmodsens_create3(pi/(nx/2)*kinfo.k.', mod_kspace, c(:,:, slcs_actual, :,:), arg.fm(:,:, slcs_actual), tts, 'mask', arg.mask(:,:, slcs_actual, :), 'nufft', arg.nufft); % "DFT and modulation matrix" times sensitivity matrix.
    
    %%% Regularization.
    if isequal(arg.beta, arg.beta(1)*ones(size(arg.beta)))
        RR = Reg1(arg.mask(:,:, slcs_actual, :), 'type_penal', 'def', 'type_diff', 'mex', 'offsets', penalty_offsets([], [nx ny]), 'order', 1, 'beta', arg.beta(1), 'distance_power', 1); % Using the same beta for each simul slice.
        R = RR.C;
        %R = reg_create3(nx, ny, nslcs, arg.nmaps, arg.beta, arg.mask(:,:, slcs_actual, :)); % For soft-SENSE.
    else
        R = reg_create1(nx, ny, nslcs, arg.beta, arg.mask(:,:, slcs_actual, :));
    end
    if dolcurves
        R1 = reg_create1(nx, ny, nslcs, 1, arg.mask(:,:, slcs_actual, :)); % Regularization fatrix for beta == 1. For L-curve computation.
    end
    
    %%% Reconstruct slices.
    % Initial solution.
    init = arg.xinit(:,:, slcs_actual);
    for map = 1:arg.nmaps; init(:,:,:, map) = init(:,:,:, 1); end
    init = init(arg.mask(:,:, slcs_actual, :));
    
    % irt toolbox CG implementation.
    if 0%save_iters
        %[xall cginfo] = qpwls_pcg(init, AC, 1, reshape(ksp_dat(:, :, acq), [], 1), 0, R, 1, arg.niter, true(size(init)), 0);
        [xall cginfo] = qpwls_pcg1(init, AC, 1, reshape(ksp_dat(:, :, acq), [], 1), R, 'niter', arg.niter, 'isave', [1:arg.niter], 'stop_diff_tol', arg.stop_diff_tol);
        x = embed(xall(:, end), reshape(arg.mask(:,:, slcs_actual, :), [], 1)); % Just want result from the last iteration.
    else
        [x cginfo] = qpwls_pcg1(init, AC, 1, reshape(ksp_dat(:, :, acq), [], 1), R, 'niter', arg.niter, 'stop_diff_tol', arg.stop_diff_tol);
        x = embed(x, reshape(arg.mask(:,:, slcs_actual, :), [], 1));
    end
    prt ['Num iters done: ' num2str(size(cginfo, 1))]
    
    % Alan's CG implementation.
    %x = conj_grad(AC, reshape(ksp_dat(:, :, acq), [], 1), init, arg.niter);
    
    imarev(:,:, slcs_actual, :) = reshape(x, nx, ny, nslcs, arg.nmaps); % Save reconstructed volumes for all acqs.
    
    if plot_extra == 2 | plot_extra == 1
        for iter = 2:2:size(cginfo, 1)
            % Plot reconstructed slices.
            xtmp = embed(xall(:, iter), arg.mask(:,:, slcs_actual, :));
            figure;
            %imagesc(tile(abs(xtmp), [nslcs 1])', clims); axis image, axis xy; cbar; colormap(gray)
            imagesc(tile(abs(xtmp), [nslcs 1])'); axis image, axis xy; cbar; colormap(gray)
            title(['Reconstructed Slices: ' num2str(slcs_actual) ', Iter ' num2str(iter)])
            pause
        end
    end
    if plot_extra == 2
        % Plot original (conventionally acquired) slices.
        load([fn_sensmap(1:end - 8) '.mat']);
        figure;
        imagesc(tile(abs(imar(:,:, slcs_actual)), [nslcs 1])'); axis image, axis xy; cbar; colormap(gray)
        title(['"Original" (conventional) Slices: ' num2str(slcs_actual)])
    end
    
    %keyboard
    
    %{
    % Calculate residual norm and regularized solution norm for the L-curve method.
    moderr = norm(AC*x(arg.mask(:,:, slcs_actual, :)) - reshape(ksp_dat(:, :, acq), [], 1), 2); % Model error.
    prerr = norm(R1*x(arg.mask(:,:, slcs_actual, :)), 2); % Prior error.
    lcur(ii, :) = [moderr prerr];
    end % for ii = 1:length(betas)
    %}
    
end % for acq = 1:nacqs
