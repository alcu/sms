%%% Reconstruction of non-SMS data.
%%%
%%% Alan Chu

%% Parameters.
do_spir = 0; % Do spiral or concentric ring k-space trajectory.
ncoils = 32;
nslcs = 39; % Overall number of slices.
nx = 64; % Object domain.
ny = nx;
ntotal_tps = 91; % Total number of time points in the run.
tps_to_do = 1;%[1:ntotal_tps]; % Recon these fMRI time points (time frame volumes).
slctd = 22;%[1:nslcs]; % Slices to transform. Use [1:nslcs] for all slices.
coiltd = [1:ncoils]; % Coil to do. Use [1:ncoils] for all coils.
lrshift = 0; % FOV shift (offset) in image domain.
tbshift = 0;
iso_slc = 20; % For readout Gz demodulation only.
opslthick = 0.3; % (cm) For readout Gz demodulation only.
use_fmap = 1; % Use field map or not.
compute_fmap = 0; % Will automatically save the field map! Compute field map, or load pre-computed fmap.
demod = 0; % Demodulate readout z-gradient from k-space data.
do_plots = 0;
do_interp = 0; % Interpolate concentric ring k-space data to see what effect interpolation has.
ninterp = 208; % Number of points to interpolate to for each segment (ring).
do_compr = 1; % Do coil compression or not.
compute_compr = 1; % Will automatically save matrices! Compute compression matrices, or load pre-computed matrices.
nvcoils = 5; % Number of virtual coils to use for coil compression.
tpcompr_to_do = 1;%[round(ntotal_tps/2) - 5:round(ntotal_tps/2) + 4]; % Time points to use for coil compression computation.

%% Experiment data locations.
% ---
fn_prefix =           '../data/ss/vol_e808_10_24_114_%04u';
fn_info =             '../data/ss/info.mat';
fn_fmap =             '../data/ss/fmfile.mat';
fn_comp_fmap_nonsms = '../data/ss/fmap_niter15_sh(0,0).mat'; % Pre-computed fmap.
%save_dir =            '../data/save_dir/';
%fn_compr =            [save_dir '/compr.mat']; % Pre-computed compression matrices.
% ---


%% Correct NFS links, if needed.
fn_list = who('fn_*', 'save_dir');
for ii = 1:length(fn_list)
    eval([fn_list{ii} ' = correct_nfs_link_path(' fn_list{ii} ');'])
end

%% Load scan info.
load(fn_info)
if do_spir
    %% Load spiral trajectory.
    kinfo.k = kinfo.k(end:-1:1); % Make trajectory go from out to in.
    kinfo.kdens = kinfo.kdens(end:-1:1);
    kinfo.kmod = kinfo.kmod(end:-1:1);
    kinfo.t = flipvec(scaninfo.te*1e-3 - kinfo.t);
    scaninfo.revflg = 3; % Use tts = kinfo.t in rec_nufft(). Also needed for rec_fmap(), since it uses rec_nufft().
    
    %% Shift trajectory for incorrect daqdel.
    del = -1;
    %del = 0;
    if del < 0
        kinfo.k = [ones([1 -del])*kinfo.k(1) kinfo.k(1:end + del)];
        kinfo.kdens = [ones([1 -del])*kinfo.kdens(1) kinfo.kdens(1:end + del)];
        kinfo.kmod = [ones([1 -del])*kinfo.kmod(1) kinfo.kmod(1:end + del)];
    elseif del > 0
        kinfo.k = [kinfo.k(del+1:end) kinfo.k(end)*ones([1 del])];
        kinfo.kdens = [kinfo.kdens(del+1:end) kinfo.kdens(end)*ones([1 del])];
        kinfo.kmod = [kinfo.kmod(del+1:end) kinfo.kmod(end)*ones([1 del])];
    end
else
    %% Load concentric rings trajectory.
    load_concen_rings_traj % Script.
end

%% FOV shift.
xsh = -1*lrshift/nx; ysh = -1*tbshift/ny;
if do_spir
    kshmod = col(exp(i*2*pi*(xsh*real(kinfo.k) + ysh*imag(kinfo.k))));
else
    kshmod = col(exp(i*2*pi*(xsh*nx/2*kx_circ/fovk_unfudged + ysh*ny/2*ky_circ/fovk_unfudged)));
end

%% Field map.
if use_fmap
    if compute_fmap
        load(fn_fmap); % Load field map data from conventional scan. Vars: ksp_dat0, ksp_dat1
        if do_spir
            kd0 = ksp_dat0;
            kd1 = ksp_dat1;
        else
            kd0 = ksp_dat0(idx_circ, :, :); % Goes from out-to-in. [ndat, ncoils, nslcs]
            kd1 = ksp_dat1(idx_circ, :, :); % Goes from out-to-in. [ndat, ncoils, nslcs]
        end
        kd0 = fov_shift(kd0, kshmod);
        kd1 = fov_shift(kd1, kshmod);
        kd0 = kd0(:, coiltd, slctd);
        kd1 = kd1(:, coiltd, slctd);
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
        kdcomp = zeros(length(idx_circ), length(coiltd), length(slctd), length(tpcompr_to_do));
        for tp_idx = 1:length(tpcompr_to_do)
            tp = tpcompr_to_do(tp_idx);
            load(sprintf(fn_prefix, tp));
            kd = ksp_dat(idx_circ, :, :); % SMS k-space data: kd goes from out-to-in. [ndat, ncoils, nacqs]
            kd = fov_shift(kd, kshmod);
            kd = kd(:, coiltd, slctd);

            kdcomp(:,:,:, tp_idx) = kd;
        end

        [kdc compr svals] = coil_compr(kdcomp, nvcoils, 0); % compr: [length(coiltd), nvcoils, length(slctd)]
        
        % Save compression matrices and singular values.
        if exist('save_dir') & isdir(save_dir)
            save([save_dir '/compr.mat'], 'compr', 'svals')
        end
    else
        load(fn_compr)
        compr = compr(coiltd, :, slctd);
    end
    
    ncoils = nvcoils;
else
    nvcoils = ncoils;
end

ntps = length(tps_to_do); % Number of fMRI time points (time frame volumes).
%xc = zeros(nx, ny, length(slctd), ntps); % Images combined across coils, for all time points.
fn_cur = mfilename;
determine_tp_vars; % Script.
start_time = clock;
tic
while true
    determine_tp; % Script.

    fprintf(['\nDoing time point ' num2str(tp) ' of ' num2str(ntps) '...\n']);

    %% Load experiment data.
    % Load k-space data and locations for 1 time point. Vars: ksp_dat
    load(sprintf(fn_prefix, tp)); % ksp_dat: [ndat, ncoils, nslcs]
    if do_spir
        kd = ksp_dat;
    else
        kd = ksp_dat(idx_circ, :, :); % k-space data: kd goes from out-to-in. [ndat, ncoils, nslcs]
    end
    kd = fov_shift(kd, kshmod);
    kd = kd(:, coiltd, slctd); % [ndat num_coiltd num_slctd]
    if do_compr
        kd = coil_compr(kd, nvcoils, 1, compr);
    end

    if exist('demod') & demod == 1
        % Compute mkss: readout Gz modulation for single-slices.
        phdiff_all = zeros(length(kinfo.t), ncoils, nslcs); % (rad) Compute for all slices, then pick ones to use based on slctd.
        count = 1;
        for slc = 1:nslcs
            slc_dist = (iso_slc - count)*opslthick; % Distance of current slice to z-gradient isocenter. (cm)
            
            % Spirals.
            %phdiff_all(:, 1, acq, slc) = gamma*1000*2*pi * cumsum(gzmod)*dt * slc_dist; % (rad) Gz modulation phase for this slice. TODO.
            
            % Concentric rings.
            phdiff_all(:, 1, slc) = 2*pi * kz_circ/100 * slc_dist; % (rad) Gz modulation phase for this slice.
            
            count = count + 1;
        end
        for coil = 2:ncoils % Make modulation the same for all coils.
            phdiff_all(:, coil, :) = phdiff_all(:, 1, :);
        end
        mkss = exp(i*phdiff_all);
        if do_compr
            mkss = mkss(:, :, slctd);
        else
            mkss = mkss(:, coiltd, slctd);
        end
    end
    if do_interp
        tts = kinfo.t;
        if exist('demod') & demod == 1
            [kd_interp_tmp kl_interp tts_interp mkss_interp] = interp_kspace(kd, kx_circ, ky_circ, tts, mkss, ninterp, length(slctd), size(kd, 2), nx); % [nx/2, ninterp, num_coiltd, num_slctd]
        else
            [kd_interp_tmp kl_interp tts_interp] = interp_kspace(kd, kx_circ, ky_circ, tts, [], ninterp, length(slctd), size(kd, 2), nx); % [nx/2, ninterp, num_coiltd, num_slctd]
        end
        kinfo_interp.k = nx/2*[col(kl_interp.').' (kx_circ(end) + i*ky_circ(end))]/fovk_unfudged;
        kinfo_interp.kmod = ones(1, nx/2*ninterp + 1);
        kinfo_interp.t = [col(tts_interp.').' tts(end)];
        g = [diff(kinfo_interp.k) 0];
        tmp_ang = unwrap(angle(kinfo_interp.k(1:5)));
        if tmp_ang(1) > tmp_ang(end) % Clockwise trajectory.
            kinfo_interp.kdens = -1*abs(g).*sin(angle(g)-angle(kinfo_interp.k));
        else % Counterclockwise trajectory.
            kinfo_interp.kdens = abs(g).*sin(angle(g)-angle(kinfo_interp.k));
        end

        kd_interp = zeros(size(kd_interp_tmp, 1)*size(kd_interp_tmp, 2) + 1, size(kd, 2), length(slctd));
        for slc = 1:length(slctd)
            for coil = 1:size(kd, 2)
                if exist('demod') & demod == 1
                    kd_interp(:, coil, slc) = [col(kd_interp_tmp(:, :, coil, slc).'); kd(end, coil, slc)] .* ...
                        exp(i * -1 * angle([col(mkss_interp(:, :, coil, slc).'); mkss(end, coil, slc)]));
                else
                    kd_interp(:, coil, slc) = [col(kd_interp_tmp(:, :, coil, slc).'); kd(end, coil, slc)];
                end
            end
        end
    else % Don't interp.
        if exist('demod') & demod == 1
            for slc = 1:length(slctd)
                for coil = 1:size(kd, 2)
                    kd(:, coil, slc) = ...
                        kd(:, coil, slc) .* ...
                        exp(i * -1 * angle(mkss(:, coil, slc)));
                end
            end
        end
    end

    %% Transform k-space data into image data.
    if do_interp
        dummy = rec_nufft(kd_interp, kinfo_interp, 'scaninfo', scaninfo, 'imsize', [nx ny], 'fm', fm, 'beta', 273, 'stop_diff_tol', 0, 'niter', 5);
    else
        dummy = rec_nufft(kd, kinfo, 'scaninfo', scaninfo, 'imsize', [nx ny], 'fm', fm, 'beta', 273, 'stop_diff_tol', 0, 'niter', 5);
    end
    x = permute(dummy, [1 2 4 3]);

    %% Combine image data from all coils into one final image.
    xc = squeeze(sqrt(sum(abs(x).^2, 3))); % x combined using sum-of-squares.

    %% Plots.
    if do_plots
        figure;
        %img(xc, 'clims', [0 400]); title(['NUFFT CG SS Recon: tp ' num2str(tp)]) % clim good for ACR phantom.
        img(xc, 'clims', [0 95]); title(['NUFFT CG SS Recon: tp ' num2str(tp)]) % clim good for human brain.
        %img(xc, 'clims', [0 11e4]); title(['Gridding and ifft2 SS Recon: tp ' num2str(tp)]) % clim good for ACR phantom.
        %img(xc, 'clims', [0 3.1e4]); title(['Gridding and ifft2 SS Recon: tp ' num2str(tp)]) % clim good for human brain.
    end

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
