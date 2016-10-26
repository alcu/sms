%%% SENSE reconstruction of simultaneous multislice data.
%%%
%%% Alan Chu

%% Parameters.
do_spir = 0; % Do spiral or concentric ring k-space trajectory.
ncoils = 32;
nslcs = 3; % Number of simultaneously acquired slices.
nslcs_overall = 39; % Overall number of slices (after separation).
nacqs = nslcs_overall/nslcs; % Number of data acquisitions for this volume.
nx = 64; % Object domain.
ny = nx;
ntotal_tps = 273; % Total number of time points in the run.
tps_to_do = 1;%[1:ntotal_tps]; % Recon these fMRI time points (time frame volumes).
acqtd = 9;%[1:nacqs]; % Acqs to recon. Use [1:nacqs] for all slices.
slctd = zeros(1, nslcs*length(acqtd));
for slc = 1:nslcs
    slctd([1:length(acqtd)] + (slc - 1)*length(acqtd)) = acqtd + (slc - 1)*nacqs;
end
lrshift = 0; % FOV shift (offset) in image domain.
tbshift = 0;
iso_slc = 20; % For readout Gz modulation. "Slice number" of z-grad isocenter location, from 1 to nslcs_overall. e.g. If iso btw slcs 2 and 3, use 2.5.
opslthick = 0.3; % (cm) For readout Gz modulation.
use_fmap = 1; % Use field map or not.
compute_fmap = 0; % Will automatically save the field map! Compute field map, or load pre-computed fmap.
compute_sensmap = 0; % Compute sensitivity map, or load pre-computed fmap.
compute_sensmap_from_ksp = 0;
doconv = 0; % Do SENSE recon for conventional, single-slice data.
do_compr = 1; % Do coil compression or not.
compute_compr = 0; % Will automatically save matrices! Compute compression matrices, or load pre-computed matrices.
nvcoils = 14; % Number of virtual coils to use for coil compression.
tpcompr_to_do = 1;%[round(ntotal_tps/2) - 5:round(ntotal_tps/2) + 4]; % Time points to use for coil compression computation.

%% Experiment data locations.
% ---
% Axial slices, 32 channels, concentric rings:
fn_prefix =           '../data/sms/vol_e808_10_24_114_%04u';
fn_info =             '../data/sms/info.mat';
fn_fmap =             '../data/sms/calib/fmfile.mat';
fn_comp_fmap_nonsms = '../data/sms/fmap_niter15_sh(0,0).mat'; % Pre-computed fmap.
fn_nonsms =           '../data/sms/calib/vol_e808_10_24_114_0021.mat';
fn_comp_nonsms =      '../data/sms/recon_x_fmap_niter15_sh(0,0)_coil14.mat'; % Pre-computed nonsms volume.
fn_comp_sensmap =     '../data/sms/sens_ESPIRiT_c_sh(0,0)_coil14.mat'; % Vars: c
%save_dir =            '../data/save_dir/'; % Directory to save vars. Don't forget trailing slash.
fn_compr =            '../data/sms/compr_coil14.mat'; % Pre-computed compression matrices.
% ---


%% Correct NFS links, if needed.
fn_list = who('fn_*', 'save_dir');
for ii = 1:length(fn_list)
    eval([fn_list{ii} ' = correct_nfs_link_path(' fn_list{ii} ');'])
end

%% Load scan info.
load(fn_info)
% For 3D-SMS recon.
fovz = nslcs_overall*opslthick/100; % (m)
resz = fovz/nslcs; % (m)
fovkz_unfudged = 1/resz/2; % Maximum extent of k-space. (cycles/m)
if do_spir
    %% Load spiral trajectory.
    kinfo.k = kinfo.k(end:-1:1); % Make trajectory go from out to in.
    kinfo.kdens = kinfo.kdens(end:-1:1);
    kinfo.kmod = kinfo.kmod(end:-1:1);
    kinfo.t = flipvec(scaninfo.te*1e-3 - kinfo.t);
    scaninfo.revflg = 3; % Use tts = kinfo.t in rec_nufft(). Also needed for rec_fmap(), since it uses rec_nufft().
    idx_circ = [1:length(kinfo.k)];
    load('spirals_data/gzmod.mat') % Vars: gzmod (Gauss/cm)
    
    %% Shift trajectory for incorrect daqdel.
    del = -1;
    %del = 0;
    if del < 0
        kinfo.k = [ones([1 -del])*kinfo.k(1) kinfo.k(1:end + del)];
        kinfo.kdens = [ones([1 -del])*kinfo.kdens(1) kinfo.kdens(1:end + del)];
        kinfo.kmod = [ones([1 -del])*kinfo.kmod(1) kinfo.kmod(1:end + del)];
        
        gzmod = [zeros(-del, 1); gzmod(1:end + del)];
    elseif del > 0
        kinfo.k = [kinfo.k(del+1:end) kinfo.k(end)*ones([1 del])];
        kinfo.kdens = [kinfo.kdens(del+1:end) kinfo.kdens(end)*ones([1 del])];
        kinfo.kmod = [kinfo.kmod(del+1:end) kinfo.kmod(end)*ones([1 del])];
        
        gzmod = [gzmod(del + 1:end); ones(del, 1)*gzmod(end)];
    end
    
    gamma = 4.2576; % (kHz/Gauss)
    dt = 4e-6; % Grad sample size. (s)
    
    % For 3D-SMS recon.
    kz_act = gamma*1000*100 * cumsum(gzmod)*dt; % (cycles/m)
    kz = nslcs/2*kz_act/fovkz_unfudged;
else
    %% Load concentric rings trajectory.
    load_concen_rings_traj % Script.
    
    % For 3D-SMS recon.
    kz = nslcs/2*kz_circ/fovkz_unfudged;
end

%% FOV shift.
xsh = -1*lrshift/nx; ysh = -1*tbshift/ny;
if do_spir
    kshmod = col(exp(i*2*pi*(xsh*real(kinfo.k) + ysh*imag(kinfo.k))));
else
    kshmod = col(exp(i*2*pi*(xsh*nx/2*kx_circ/fovk_unfudged + ysh*ny/2*ky_circ/fovk_unfudged)));
end

%%% Modulations.
phdiff_all = zeros(length(kinfo.k), ncoils, nacqs, nslcs); % (rad) Compute for all slices, then pick ones to use based on slctd.
if ~exist('doconv') | doconv == 0
    count = 1;
    for slc = 1:nslcs
        for acq = 1:nacqs
            slc_dist = (iso_slc - count)*opslthick; % Distance of current slice to z-gradient isocenter. (cm)
            
            if do_spir
                % Spirals.
                phdiff_all(:, 1, acq, slc) = gamma*1000*2*pi * cumsum(gzmod)*dt * slc_dist; % (rad) Gz modulation phase for this slice. TODO.
            else
                % Concentric rings.
                phdiff_all(:, 1, acq, slc) = 2*pi * kz_circ/100 * slc_dist; % (rad) Gz modulation phase for this slice.
            end
            
            count = count + 1;
        end
    end
    for coil = 2:ncoils % Make modulation the same for all coils.
        phdiff_all(:, coil, :, :) = phdiff_all(:, 1, :, :);
    end
end
phdiff_all = phdiff_all(:,:, acqtd,:);

%%% Field map.
if use_fmap
    if compute_fmap
        %%% Concentric rings kx-ky trajectory.
        load(fn_fmap); % Load field map data from conventional scan. Vars: ksp_dat0, ksp_dat1
        kd0 = ksp_dat0(idx_circ, :, :); % Goes from out-to-in. [ndat, ncoils, nslcs_overall]
        kd1 = ksp_dat1(idx_circ, :, :); % Goes from out-to-in. [ndat, ncoils, nslcs_overall]
        kd0 = fov_shift(kd0, kshmod);
        kd1 = fov_shift(kd1, kshmod);
        kd0 = kd0(:,:, slctd);
        kd1 = kd1(:,:, slctd);

        tic
        fm = rec_fmap(kd0, kd1, kinfo, 'imsize', [nx ny], 'scaninfo', scaninfo, 'stop_diff_tol', 0, 'niter', 15);
        time_fmap = toc

        %keyboard % Save fm for later use.
        if exist('fn_comp_fmap_nonsms')
            save(fn_comp_fmap_nonsms, 'fm')
        end
    else
        load(fn_comp_fmap_nonsms); % Vars: fm
        fm = fm(:,:, slctd);
    end
else
    fm = zeros(nx, ny, length(slctd)); % No field map.
end

%% Coil compression.
if do_compr
    if compute_compr
        kdcomp = zeros(length(idx_circ), ncoils, length(acqtd), length(tpcompr_to_do));
        for tp_idx = 1:length(tpcompr_to_do)
            tp = tpcompr_to_do(tp_idx);
            load(sprintf(fn_prefix, tp));
            kd = ksp_dat(idx_circ, :, :); % SMS k-space data: kd goes from out-to-in. [ndat, ncoils, nacqs]
            kd = fov_shift(kd, kshmod);
            kd = kd(:,:, acqtd);

            kdcomp(:,:,:, tp_idx) = kd;
        end

        [kdc compr svals] = coil_compr(kdcomp, nvcoils, 0); % compr: [ncoils, nvcoils, length(acqtd)]
        
        % Save compression matrices and singular values.
        if exist('save_dir') & isdir(save_dir)
            save([save_dir '/compr.mat'], 'compr', 'svals')
        end
    else
        load(fn_compr)
        compr = compr(:, :, acqtd);
    end

    ncoils = nvcoils;
else
    nvcoils = ncoils;
end

%% Coil sensitivities.
if compute_sensmap
    if compute_sensmap_from_ksp
        %% Reconstruct 1 non-simultaneous slice volume for mask and sens map construction.
        % Load k-space data and locations for the first time point. Vars: ksp_dat
        %load(fn_nonsms);
        if ~compute_fmap | ~use_fmap
            load(fn_fmap); % Load field map data from conventional scan. Vars: ksp_dat0, ksp_dat1
        end
        ksp_dat = ksp_dat1; % Use second, non-advanced/delayed time frame of field map acqs, since it doesn't have z-gradient modulation.

        %%% Concentric rings kx-ky trajectory.
        kd = ksp_dat(idx_circ, :, :);
        kd = fov_shift(kd, kshmod);
        kd = kd(:,:, slctd);
        if do_compr
            kd = coil_compr(kd, nvcoils, 1, repmat(compr, [1 1 nslcs]));
        end

        imarall = rec_nufft(kd, kinfo, 'scaninfo', scaninfo, 'imsize', [nx ny], 'fm', fm, 'beta', 273, 'stop_diff_tol', 0, 'niter', 15);
        imarall = permute(imarall, [1 2 4 3]);
    else
        load(fn_comp_nonsms); % Vars: x
        imarall = x;
        imarall = x(:,:,:, slctd);
    end

    betasens_log2 = 2;
    imar_body = sqrt(sum(abs(imarall).^2, 3)); % Sim body coil image.
    c = zeros(nx, ny, length(slctd), ncoils);
    for slc = 1:length(slctd)
        fprintf(['\nDoing sens maps for slice ' num2str(slc) ' of ' num2str(length(slctd)) '...\n']);
        if do_compr
            coil_for_ph = 1; % Good for SVD coil compression. (1st coil has the most signal.)
        else
            coil_for_ph = 3; % Good for our 32-ch coil. (1st coil doesn't have much signal usually.)
        end
        c(:,:, slc, :) = mri_sensemap_denoise(squeeze(imarall(:,:,:, slc)), 'bodycoil', imar_body(:,:, slc).*exp(i*angle(imarall(:,:, coil_for_ph, slc))), 'chol', 1, 'niter', 1, 'l2b', betasens_log2, 'order', 2); % Uses phase of 3rd coil. If don't specify 'bodycoil', default is to use 1st coil, which sucks for our 32-ch coil for some reason.
        %c(:,:, slc, :) = mri_sensemap_denoise(squeeze(imarall(:,:,:, slc)), 'chol', 1, 'niter', 1, 'l2b', betasens_log2, 'order', 2);
    end
else
    load(fn_comp_sensmap);
    c = c(:,:, slctd,:, 2);
end


ntps = length(tps_to_do);
%imarev = zeros(nx, ny, nslcs_overall, ntps); % Save reconstructed volumes for all time points.
fn_cur = mfilename;
determine_tp_vars; % Script.
start_time = clock;
tic
while true
    determine_tp; % Script.

    fprintf(['\nDoing time point ' num2str(tp) ' of ' num2str(ntps) '...\n']);

    % Load k-space data and locations for 1 time point. Vars: ksp_dat
    load(sprintf(fn_prefix, tp));

    %%% Concentric rings kx-ky trajectory.
    kd = ksp_dat(idx_circ, :, :);
    kd = fov_shift(kd, kshmod);
    kd = kd(:,:, acqtd);
    if do_compr
        kd = coil_compr(kd, nvcoils, 1, compr);
    end

    xc = rec_sense(kd, kinfo, c, 'imsize', [nx ny], 'fm', fm, 'gzmod', phdiff_all, 'scaninfo', scaninfo, 'beta', [273 273 273], 'niter', 10, 'stop_diff_tol', 0);
    %xc = rec_3d(kd, kinfo, c, kz, nacqs, slctd, 'imsize', [nx ny], 'fm', fm, 'scaninfo', scaninfo, 'beta', [273 273 273], 'niter', 10, 'stop_diff_tol', 0);

    % Save some vars in case of crash.
    if exist('save_dir') & isdir(save_dir)
        save([save_dir '/recon_tp' num2str(tp, '%.3d') '.mat'], 'xc')
    end

    % Output progress.
    cur_time = clock;
    elap_time = etime(cur_time, start_time); % (s)
    fileID = fopen(['/tmp/' fn_cur '_progress.txt'], 'w');
    fprintf(fileID, 'Finished time point %d of %d.\n', tp_idx, ntps);
    fprintf(fileID, 'Elapsed time is %f hours.\n', elap_time/60/60);
    fprintf(fileID, 'Estimated time remaining is %f hours.\n', (elap_time / (tp_idx/ntps)  - elap_time)/60/60);
    fclose(fileID);
end
time_recon = toc

figure;
img(xc) % Reconstructed image.
