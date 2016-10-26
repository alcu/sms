%%% GRAPPA-like recon of non-Cartesian simultaneous multislice data.
%%%
%%% Alan Chu

%% Parameters.
ncoils = 32;
nslcs = 3; % Number of simultaneously acquired slices.
nslcs_overall = 39; % Overall number of slices (after separation).
nacqs = nslcs_overall/nslcs; % Number of data acquisitions for this volume.
nx = 64; % Object domain.
ny = nx;
ntotal_tps = 273; % Total number of time points in the run.
tps_to_do = 1;%[1:ntotal_tps]; % Recon these fMRI time points (time frame volumes).
ave_to_do = 21;%[round(ntotal_tps/2)]; % Time points over which to average GRAPPA weight computation.
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
do_split = 1; % Do split slice-GRAPPA. ave_to_do has to be only 1 time frame.
compute_kernel = 1; % Will automatically save the kernel! Compute GRAPPA kernel, or load pre-computed kernel.
do_compr = 1; % Do coil compression or not.
do_compr_ss = 1; % If doing coil compression, do compression on SS (target) only, using GRAPPA to do SMS coil compression.
compute_compr = 1; % Will automatically save matrices! Compute compression matrices, or load pre-computed matrices.
nvcoils = 5; % Number of virtual coils to use for coil compression.
ntotal_tps_ss = 21; % Total number of time points in the SS run.
if do_compr_ss
    %tpcompr_to_do = [round(ntotal_tps_ss/2) - 5:round(ntotal_tps_ss/2) + 4]; % Time points to use for coil compression computation.
    %tpcompr_to_do = [ntotal_tps_ss - 9:ntotal_tps_ss]; % For using calibration time frames for coil compression matrix computation.
    tpcompr_to_do = 21;
else
    %tpcompr_to_do = [round(ntotal_tps/2) - 5:round(ntotal_tps/2) + 4]; % Time points to use for coil compression computation.
    tpcompr_to_do = 1;
end
do_plots = 0;
do_test_trans = 0;
ncirccent = 3; % Number of rings around center k-space point to use for center weight computation.
ninterp = 208; % Number of points to interpolate to for each segment (ring).
nradsects = 1;
nangsects = 8; % Make sure ninterp is divisible by this.
%% Params for 1st (outer) radial sector.
nrings(1) = 32; % Number of rings to use for this radial sector.
xkern_coord{1} = [1 2 3]; % x kernel weight coordinates (can "skip over" certain source positions). [1, nrows]
ykern_coord{1} = [1 2 3]; % y kernel weight coordinates (can "skip over" certain source positions). [1, ncols]
kern_cent{1} = [xkern_coord{1}(round(length(xkern_coord{1})/2)) ykern_coord{1}(round(length(ykern_coord{1})/2))]; % Kernel center (coordinate for where the kernel computes the value). [1, 2]
nradovlap{1} = kern_cent{1}(1) - 1; % Number of samples that overlap beyond the radial sector boundary.
rad_trans_wts{1} = [zeros(1, nradovlap{1}) ones(1, nradovlap{1})]; % No smoothing. [1, 2*nradovlap{radsect}]
nangovlap{1} = kern_cent{1}(2) - 1; % Number of samples that overlap beyond the angular sector boundary.
ang_trans_wts{1} = [zeros(1, nangovlap{1}) ones(1, nangovlap{1})]; % No smoothing. [1, 2*nangovlap{radsect}]
% ang_trans_wts{1} = 1/(2*nangovlap{1} - (2*(kern_cent{1}(2) - 1) - 1)) * ...
%     [zeros(1, kern_cent{1}(2) - 2), ...
%      [0:1:2*nangovlap{1} - (2*(kern_cent{1}(2) - 1) - 1)], ...
%      (2*nangovlap{1} - (2*(kern_cent{1}(2) - 1) - 1))*ones(1, kern_cent{1}(2) - 2)];
%% Params for 2nd radial sector.
nrings(2) = 16; % Number of rings to use for this radial sector.
xkern_coord{2} = [1 2 3 4 5]; % x kernel weight coordinates (can "skip over" certain source positions). [1, nrows]
ykern_coord{2} = [1 3 5 7 9]; % y kernel weight coordinates (can "skip over" certain source positions). [1, ncols]
kern_cent{2} = [xkern_coord{2}(round(length(xkern_coord{2})/2)) ykern_coord{2}(round(length(ykern_coord{2})/2))]; % Kernel center (coordinate for where the kernel computes the value). [1, 2]
nradovlap{2} = kern_cent{2}(1) - 1; % Number of samples that overlap beyond the radial sector boundary.
rad_trans_wts{2} = [zeros(1, nradovlap{2}) ones(1, nradovlap{2})]; % No smoothing. [1, 2*nradovlap{radsect}]
nangovlap{2} = kern_cent{2}(2) - 1; % Number of samples that overlap beyond the angular sector boundary.
ang_trans_wts{2} = [zeros(1, nangovlap{2}) ones(1, nangovlap{2})]; % No smoothing. [1, 2*nangovlap{radsect}]
% ang_trans_wts{2} = 1/(2*nangovlap{2} - (2*(kern_cent{2}(2) - 1) - 1)) * ...
%     [zeros(1, kern_cent{2}(2) - 2), ...
%      [0:1:2*nangovlap{2} - (2*(kern_cent{2}(2) - 1) - 1)], ...
%      (2*nangovlap{2} - (2*(kern_cent{2}(2) - 1) - 1))*ones(1, kern_cent{2}(2) - 2)];
%% Params for 3rd radial sector.
nrings(3) = 6; % Number of rings to use for this radial sector.
xkern_coord{3} = [1 2 3 4 5]; % x kernel weight coordinates (can "skip over" certain source positions). [1, nrows]
ykern_coord{3} = [1 5 9 13 17]; % y kernel weight coordinates (can "skip over" certain source positions). [1, ncols]
kern_cent{3} = [xkern_coord{3}(round(length(xkern_coord{3})/2)) ykern_coord{3}(round(length(ykern_coord{3})/2))]; % Kernel center (coordinate for where the kernel computes the value). [1, 2]
nradovlap{3} = kern_cent{3}(1) - 1; % Number of samples that overlap beyond the radial sector boundary.
rad_trans_wts{3} = [zeros(1, nradovlap{3}) ones(1, nradovlap{3})]; % No smoothing. [1, 2*nradovlap{radsect}]
nangovlap{3} = kern_cent{3}(2) - 1; % Number of samples that overlap beyond the angular sector boundary.
ang_trans_wts{3} = [zeros(1, nangovlap{3}) ones(1, nangovlap{3})]; % No smoothing. [1, 2*nangovlap{radsect}]

%% Experiment data locations.
% ---
fn_prefix =           '../data/sms/vol_e808_10_24_114_%04u';
fn_nonsms_gzmod =     '../data/sms/calib/vol_e808_10_24_114_0021.mat';
fn_prefix_sms_gzmod = '../data/sms/calib/vol_summed_%04u';
fn_info =             '../data/sms/info.mat';
fn_fmap =             '../data/sms/calib/fmfile.mat';
fn_comp_fmap_nonsms = '../data/sms/fmap_niter15_sh(0,0).mat'; % Pre-computed fmap.
fn_ss_compr =         '../data/sms/calib/vol_e808_10_24_114_%04u';
%save_dir =            '../data/save_dir/'; % Directory to save vars. Don't forget trailing slash.
%fn_kernel =           [save_dir '/kernel.mat']; % Pre-computed GRAPPA kernel.
%fn_compr =            [save_dir '/compr.mat']; % Pre-computed compression matrices.
% ---


%% Error checks and exceptions.
if sum(nrings(1:nradsects)) ~= nx/2
    error('The number of rings for each radial sector is incorrect. Check value(s) of ''nrings''.')
end
if nradsects == 1
    nradovlap{1} = 0; % nradovlap{1} must be 0 if nradsects == 1 because can't do any radial overlaps.
end
if mod(ninterp, nangsects)
    error('''ninterp'' is not divisible by ''nangsects''.')
end

%% Correct NFS links, if needed.
fn_list = who('fn_*', 'save_dir');
for ii = 1:length(fn_list)
    eval([fn_list{ii} ' = correct_nfs_link_path(' fn_list{ii} ');'])
end

%% Load scan info and concentric rings trajectory.
load(fn_info)
load_concen_rings_traj % Script.

%% FOV shift.
xsh = -1*lrshift/nx; ysh = -1*tbshift/ny;
kshmod = col(exp(i*2*pi*(xsh*nx/2*kx_circ/fovk_unfudged + ysh*ny/2*ky_circ/fovk_unfudged)));

%% Field map.
if use_fmap
    if compute_fmap
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
        %figure; img(fm, 'clims', [-225 225]); title('fmap')
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
        if do_compr_ss
            kdcomp = zeros(length(idx_circ), ncoils, length(slctd), length(tpcompr_to_do));
        else
            kdcomp = zeros(length(idx_circ), ncoils, length(acqtd), length(tpcompr_to_do));
        end
        for tp_idx = 1:length(tpcompr_to_do)
            tp = tpcompr_to_do(tp_idx);
            if do_compr_ss
                load(sprintf(fn_ss_compr, tp));
            else
                load(sprintf(fn_prefix_sms_gzmod, tp));
            end
            kd = ksp_dat(idx_circ, :, :); % SMS k-space data: kd goes from out-to-in. [ndat, ncoils, nacqs]
            kd = fov_shift(kd, kshmod);
            if do_compr_ss
                kd = kd(:,:, slctd);
            else
                kd = kd(:,:, acqtd);
            end

            kdcomp(:,:,:, tp_idx) = kd;
        end

        [kdc compr svals] = coil_compr(kdcomp, nvcoils, 0); % compr: [ncoils, nvcoils, length(acqtd) OR length(slctd)]
        
        % Save compression matrices and singular values.
        if exist('save_dir') & isdir(save_dir)
            save([save_dir '/compr.mat'], 'compr', 'svals')
        end
    else
        load(fn_compr)
        if do_compr_ss
            compr = compr(:, :, slctd);
        else
            compr = compr(:, :, acqtd);
        end
    end
else
    nvcoils = ncoils;
end


if compute_kernel
    do_wts = 1;
else
    do_wts = 0;
    load(fn_kernel); % Vars: w_cent, w_sect, w_sect_asym
    for radsect = 1:nradsects
        for angsect = 1:nangsects
            w_sect{radsect, angsect} = w_sect{radsect, angsect}(:,:, acqtd);
        end
    end
    w_cent = w_cent(:,:, acqtd);
    for edg = 1:2 % Inner and outer edges.
        for angsect = 1:nangsects
            w_sect_asym{edg, angsect} = w_sect_asym{edg, angsect}(:,:, acqtd);
        end
    end
end
do_calib = 1;
do_recon = 0;
nave = length(ave_to_do);
fn_cur = mfilename;
first_tp = 1;
tic
for tp_idx = 1:nave
    tp = ave_to_do(tp_idx);
    tp_idx_ave = tp_idx;

    fprintf(['\nDoing time point ' num2str(tp) ' of ' num2str(nave) '...\n']);

    grappa_concen_rings % Script.
end
time_calib = toc
% Save kernel.
if exist('save_dir') & isdir(save_dir) & do_wts == 1
    save([save_dir '/kernel.mat'], 'w_sect', 'w_sect_asym', 'w_cent')
end


do_wts = 0;
do_calib = 0;
do_recon = 1;
ntps = length(tps_to_do);
%xc = zeros(nx, ny, nslcs_overall, ntps); % Images combined across coils, for all time points.
%xc = zeros(nx, ny, ntps); % Special case: Save just 1 coil of 1 separated slice.
%xc = zeros(nx, ny, ncoils, ntps); % Special case: Save all coils of 1 separated slice.
first_tp = 1;
determine_tp_vars; % Script.
start_time = clock;
tic
while true
    determine_tp; % Script.

    fprintf(['\nDoing time point ' num2str(tp) ' of ' num2str(ntps) '...\n']);

    grappa_concen_rings % Script.

    % Save some vars in case of crash.
    if exist('save_dir') & isdir(save_dir)
        %xc_tmp = xc(:,:,:, tp);
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
