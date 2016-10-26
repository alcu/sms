%% Load experiment data.
% Load k-space data and locations for 1 time point. Vars: ksp_dat
if (do_calib == 1) & first_tp
    % Load k-space data for 1 non-simultaneous slice volume that has readout Gz modulation. Vars: ksp_dat
    load(fn_nonsms_gzmod);
    kdss = ksp_dat(idx_circ, :, :); % Single-slice k-space data: kdss goes from out-to-in. [ndat, ncoils, nslcs_overall]
    kdss = fov_shift(kdss, kshmod);
    kdss = kdss(:,:, slctd);

    tts = kinfo.t; % tts goes from out-to-in.
    % Compute mkss: readout Gz modulation for single-slices.
    phdiff_all = zeros(length(tts), nvcoils, nslcs_overall); % (rad) Compute for all slices, then pick ones to use based on slctd.
    predemod = zeros(length(tts), ncoils, nslcs_overall); % (rad) Compute for all slices, then pick ones to use based on slctd.
    count = 1;
    for slc = 1:nslcs_overall
        pre_slc_dist = (iso_slc - count)*opslthick; % Distance of current slice to z-gradient isocenter. (cm)
        slc_dist = (iso_slc - (round(nacqs/2) + ...
                               (ceil(slc/nacqs) - 1)*nacqs)) * ...
            opslthick; % Equivalent distance of current slice to isocenter after pre-demodulation.

        % Spirals.
        %phdiff_all(:, 1, acq, slc) = gamma*1000*2*pi * cumsum(gzmod)*dt * slc_dist; % (rad) Gz modulation phase for this slice. TODO.
        
        % Concentric rings.
        phdiff_all(:, 1, slc) = 2*pi * kz_circ/100 * slc_dist; % (rad) Gz modulation phase for this slice.
        predemod(:, 1, slc) = 2*pi * kz_circ/100 * pre_slc_dist - phdiff_all(:, 1, slc).'; % (rad) To make the modulation like isocenter modulation.
        
        count = count + 1;
    end
    for coil = 2:nvcoils % Make modulation the same for all coils.
        phdiff_all(:, coil, :) = phdiff_all(:, 1, :);
    end
    for coil = 2:ncoils % Make modulation the same for all coils.
        predemod(:, coil, :) = predemod(:, 1, :);
    end
    mkss = exp(i*phdiff_all);
    mkss = mkss(:,:, slctd);
    pre_mkss = exp(i*predemod);
    pre_mkss = pre_mkss(:,:, slctd);

    % Pre-demodulation.
    for coil = 1:ncoils
        for slc = 1:length(slctd)
            kdss(:, coil, slc) = kdss(:, coil, slc) .* ...
                exp(i * -1 * angle(pre_mkss(:, coil, slc)));
        end
    end

    if do_compr
        if do_compr_ss
            if do_split
                kdss_orig = kdss;
            end
            kdss = coil_compr(kdss, nvcoils, 1, compr);
        else
            kdss = coil_compr(kdss, nvcoils, 1, repmat(compr, [1 1 nslcs]));
        end
    end
end
if do_recon | (do_wts & ~do_split)
    if do_calib == 1
        load(sprintf(fn_prefix_sms_gzmod, tp)); % Load calibration SMS volume (with TR same as SS scan, and FA similar too).
    elseif (do_recon == 1) & (do_calib == 0)
        load(sprintf(fn_prefix, tp));
    end
    kd = ksp_dat(idx_circ, :, :); % SMS k-space data: kd goes from out-to-in. [ndat, ncoils, nacqs]
    kd = fov_shift(kd, kshmod);
    kd = kd(:,:, acqtd);

    % Pre-demodulation.
    for coil = 1:ncoils
        for acq = 1:length(acqtd)
            kd(:, coil, acq) = kd(:, coil, acq) .* ...
                exp(i * -1 * angle(pre_mkss(:, coil, acq)));
        end
    end

    if do_compr
        if ~do_compr_ss
            kd = coil_compr(kd, nvcoils, 1, compr);
        end
    end
end

%% Interpolate k-space data (both SMS and single-slice) to constant angular velocity trajectory.
if do_recon | (do_wts & ~do_split)
    fprintf(['\n\nInterpolating SMS k-space data...'])
    [kd_interp] = interp_kspace(kd, kx_circ, ky_circ, [], [], ninterp, length(acqtd), size(kd, 2), nx); % kd_interp: [nx/2, ninterp, ncoils, length(acqtd)]. kl_interp, tts_interp: [nx/2, ninterp].
end
if (do_calib == 1) & first_tp
    fprintf(['\n\nInterpolating single-slice k-space data...'])
    [kdss_interp kl_interp tts_interp mkss_interp] = interp_kspace(kdss, kx_circ, ky_circ, tts, mkss, ninterp, length(slctd), size(kdss, 2), nx); % [nx/2, ninterp, nvcoils, length(slctd)]
    if do_wts & do_compr & do_compr_ss & do_split
        [kdss_orig_interp] = interp_kspace(kdss_orig, kx_circ, ky_circ, [], [], ninterp, length(slctd), size(kdss_orig, 2), nx); % [nx/2, ninterp, ncoils, length(slctd)]
    end
end

%% Separate interpolated, constant ang velocity k-space data (both SMS and SS) into sectors.
%kd_sect = zeros(nradsects, nangsects, nx/2/nradsects, ninterp/nangsects, ncoils, length(acqtd));
if do_recon | (do_wts & ~do_split)
    kd_sect = cell(nradsects, nangsects); % Need cell b/c overlapping sectors will have different sizes.
end
if (do_calib == 1) & first_tp
    %kdss_sect = zeros(nradsects, nangsects, nx/2/nradsects, ninterp/nangsects, nvcoils, length(slctd));
    kdss_sect = cell(nradsects, nangsects);
    if do_wts & do_compr & do_compr_ss & do_split
        kdss_orig_sect = cell(nradsects, nangsects);
    end
end
for radsect = 1:nradsects
    for angsect = 1:nangsects
        if radsect == 1
            rad_range = [1 : nrings(radsect) + nradovlap{radsect}];
        elseif radsect == nradsects
            rad_range = [1 - nradovlap{radsect} : nrings(radsect)] + sum(nrings(1:radsect - 1));
        else % Don't need to worry about boundaries.
            rad_range = [1 - nradovlap{radsect} : nrings(radsect) + nradovlap{radsect}] + sum(nrings(1:radsect - 1));
        end
        if angsect == 1
            ang_range = [[ninterp - nangovlap{radsect} + 1 : ninterp] [1 : ninterp/nangsects + nangovlap{radsect}]];
        elseif angsect == nangsects
            ang_range = [[[1 - nangovlap{radsect} : ninterp/nangsects] + (angsect - 1)*ninterp/nangsects] [1:nangovlap{radsect}]];
        else % Don't need to worry about boundaries.
            ang_range = [1 - nangovlap{radsect} : ninterp/nangsects + nangovlap{radsect}] + (angsect - 1)*ninterp/nangsects;
        end
        if do_recon | (do_wts & ~do_split)
            kd_sect{radsect, angsect} = kd_interp(rad_range, ang_range, :, :);
        end
        if (do_calib == 1) & first_tp
            kdss_sect{radsect, angsect} = kdss_interp(rad_range, ang_range, :, :);
            if do_wts & do_compr & do_compr_ss & do_split
                kdss_orig_sect{radsect, angsect} = kdss_orig_interp(rad_range, ang_range, :, :);
            end
        end
    end
end

%% Compute GRAPPA weights and reconstruct single-slice k-space data.
%xk_sect = zeros(size(kdss_sect)); % [nradsects, nangsects, nx/2/nradsects, ninterp/nangsects, nvcoils, length(slctd)]
xk_sect = cell(nradsects, nangsects);
if first_tp
    pos_for_w_cent = [0:7]*round(ninterp)/8 + 1; % Use only these positions in each ring for computing the center k-space point.
end
if (do_wts == 1) & first_tp
    if ~do_split
        kd_sect_firstnave = cell(nradsects, nangsects); % Each of these cells will contain a cell of nave kd_sect's.
    end
    kd_interp_cent_firstnave = cell(1, nave); % Each of these cells contain a reshaped kd_interp(end - ncirccent + 1:end, pos_for_w_cent, :, :).
    if do_compr & ~do_compr_ss % Ugly hack.
        ncoils_src = nvcoils;
    else
        ncoils_src = ncoils;
    end
    for radsect = 1:nradsects
        for angsect = 1:nangsects
            nsrc = length(xkern_coord{radsect})*length(ykern_coord{radsect}); % Make sure this is the same as in grappa_weights() and grappa_recon().
            w_sect{radsect, angsect} = zeros(ncoils_src*nsrc, nvcoils*nslcs, length(acqtd));
            if radsect == 1
                w_sect_asym{1, angsect} = zeros(ncoils_src*nsrc, nvcoils*nslcs, length(acqtd)); % Asymmetric kernel weights (for outer and inner edges).
            end
            if radsect == nradsects
                w_sect_asym{2, angsect} = zeros(ncoils_src*nsrc, nvcoils*nslcs, length(acqtd)); % Asymmetric kernel weights (for outer and inner edges).
            end
            if ~do_split
                kd_sect_firstnave{radsect, angsect} = cell(1, nave);
            end
        end
    end
end
for radsect = 1:nradsects
    for angsect = 1:nangsects
        xtr = [kern_cent{radsect}(1):size(kdss_sect{radsect, angsect}, 1) - kern_cent{radsect}(1) + 1]; % x training range. [1, num_x_training_pts]
        %xtr = [18:size(kdss_sect{radsect, angsect}, 1) - kern_cent{radsect}(1) + 1]; % x training range. [1, num_x_training_pts]
        ytr = [kern_cent{radsect}(2):size(kdss_sect{radsect, angsect}, 2) - kern_cent{radsect}(2) + 1]; % y training range. [1, num_y_training_pts]
        
        xrec = [kern_cent{radsect}(1):size(kdss_sect{radsect, angsect}, 1) - kern_cent{radsect}(1) + 1]; % x recon range. [1, num_x_recon_pts]
        yrec = [kern_cent{radsect}(2):size(kdss_sect{radsect, angsect}, 2) - kern_cent{radsect}(2) + 1]; % y recon range. [1, num_y_recon_pts]
        xrec_outer = xrec(1:kern_cent{radsect}(1) - 1); % x recon range for asymmetric kernel for outer radial sector. Only used for recon.
        xrec_inner = xrec(end - kern_cent{radsect}(1) + 2:end); % x recon range for asymmetric kernel for inner radial sector. Only used for recon.
        if do_wts == 1
            fprintf(['\nSaving interpolated SMS k-space of radsect ' num2str(radsect) ' of ' num2str(nradsects) ', angsect ' num2str(angsect) ' of ' num2str(nangsects) '...\n']);
            if ~do_split
                kd_sect_firstnave{radsect, angsect}{tp_idx_ave} = kd_sect{radsect, angsect};
            end
            if tp == ave_to_do(end)
                fprintf(['\nComputing GRAPPA weights of radsect ' num2str(radsect) ' of ' num2str(nradsects) ', angsect ' num2str(angsect) ' of ' num2str(nangsects) '...\n']);
                if do_split
                    if do_compr & do_compr_ss
                        w_sect{radsect, angsect} = grappa_weights_split(nslcs, kdss_orig_sect{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr);
                    else
                        w_sect{radsect, angsect} = grappa_weights_split(nslcs, kdss_sect{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr);
                    end
                else
                    w_sect{radsect, angsect} = grappa_weights_ave(kd_sect_firstnave{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr);
                end
                
                if radsect == 1 % Compute asymmetric kernel weights for outer and inner radial sectors.
                    if do_split
                        if do_compr & do_compr_ss
                            w_sect_asym{1, angsect} = grappa_weights_split(nslcs, kdss_orig_sect{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr, -1*(kern_cent{radsect}(1) - 1));
                        else
                            w_sect_asym{1, angsect} = grappa_weights_split(nslcs, kdss_sect{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr, -1*(kern_cent{radsect}(1) - 1));
                        end
                    else
                        w_sect_asym{1, angsect} = grappa_weights_ave(kd_sect_firstnave{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr, -1*(kern_cent{radsect}(1) - 1));
                    end
                end
                if radsect == nradsects
                    if do_split
                        if do_compr & do_compr_ss
                            w_sect_asym{2, angsect} = grappa_weights_split(nslcs, kdss_orig_sect{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr, kern_cent{radsect}(1) - 1);
                        else
                            w_sect_asym{2, angsect} = grappa_weights_split(nslcs, kdss_sect{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr, kern_cent{radsect}(1) - 1);
                        end
                    else
                        w_sect_asym{2, angsect} = grappa_weights_ave(kd_sect_firstnave{radsect, angsect}, kdss_sect{radsect, angsect}, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xtr, ytr, kern_cent{radsect}(1) - 1);
                    end
                end
            end
        end
        if do_recon == 1
            fprintf(['\nComputing GRAPPA recon of radsect ' num2str(radsect) ' of ' num2str(nradsects) ', angsect ' num2str(angsect) ' of ' num2str(nangsects) '...\n']);
            xk_sect{radsect, angsect} = grappa_recon(kd_sect{radsect, angsect}, w_sect{radsect, angsect}, nslcs, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xrec, yrec);
            
            % Expand matrix with zeros to match size of kd_sect.
            xk_tmp = xk_sect{radsect, angsect};
            xk_sect{radsect, angsect} = zeros(size(kd_sect{radsect, angsect}, 1), size(kd_sect{radsect, angsect}, 2), nvcoils, length(slctd));
            xk_sect{radsect, angsect}(xrec, yrec, :, :) = xk_tmp;
            
            if radsect == 1 % Recon using asymmetric kernel weights for outer and inner radial sectors.
                xk_tmp = grappa_recon(kd_sect{radsect, angsect}, w_sect_asym{1, angsect}, nslcs, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xrec_outer, yrec);
                
                xk_sect{radsect, angsect}(1:kern_cent{radsect}(1) - 1, yrec, :, :) = xk_tmp; % Replace outermost ring(s).
            end
            if radsect == nradsects
                xk_tmp = grappa_recon(kd_sect{radsect, angsect}, w_sect_asym{2, angsect}, nslcs, xkern_coord{radsect}, ykern_coord{radsect}, kern_cent{radsect}, xrec_inner, yrec);
                
                xk_sect{radsect, angsect}(end - kern_cent{radsect}(1) + 2:end, yrec, :, :) = xk_tmp; % Replace innermost ring(s).
            end
        end
    end
end
% Special kernel for center k-space sample.
if do_wts == 1
    if do_split
        if do_compr & do_compr_ss
            kd_interp_cent = reshape(kdss_orig_interp(end - ncirccent + 1:end, pos_for_w_cent, :, :), [], size(kdss_orig_interp, 3), length(slctd));
        else
            kd_interp_cent = reshape(kdss_interp(end - ncirccent + 1:end, pos_for_w_cent, :, :), [], size(kdss_interp, 3), length(slctd));
        end
    else
        kd_interp_cent_firstnave{tp_idx_ave} = reshape(kd_interp(end - ncirccent + 1:end, pos_for_w_cent, :, :), [], size(kd_interp, 3), length(acqtd));
    end
    if tp == ave_to_do(end)
        if do_split
            w_cent = grappa_weights_center_split(nslcs, kd_interp_cent, ...
                                                 squeeze(kdss(end, :, :)));
        else
            w_cent = grappa_weights_center_ave(kd_interp_cent_firstnave, ...
                                               squeeze(kdss(end, :, :)));
        end
    end
end
if do_recon == 1
    xk_cent = grappa_recon_center(reshape(kd_interp(end - ncirccent + 1:end, pos_for_w_cent, :, :), [], size(kd_interp, 3), length(acqtd)), ...
                                  w_cent, nslcs);
end

if do_recon == 1
    %%% Combine all sectors of single-slice k-space data.
    %% First combine and smooth borders of angular sectors.
    xktmp = cell(1, nradsects); % Contains combined angular sectors with smoothed borders.
    for radsect = 1:nradsects
        for angsect = 1:nangsects
            tmp = xk_sect{radsect, angsect}(:, nangovlap{radsect} + 1 : end - nangovlap{radsect}, :,:);

            if angsect == 1
                tmp1 = xk_sect{radsect, nangsects}; % Last sector (will use as previous sector).
                tmp2 = xk_sect{radsect, angsect + 1}; % Next sector.
            elseif angsect == nangsects
                tmp1 = xk_sect{radsect, angsect - 1}; % Previous sector.
                tmp2 = xk_sect{radsect, 1}; % First sector (will use as next sector).
            else % Don't need to worry about boundaries.
                tmp1 = xk_sect{radsect, angsect - 1}; % Previous sector.
                tmp2 = xk_sect{radsect, angsect + 1}; % Next sector.
            end

            for coil = 1:nvcoils
                for slc = 1:length(slctd)
                    % Adjust beginning border.
                    tmp(:, 1:nangovlap{radsect}, coil, slc) = ...
                        repmat(ang_trans_wts{radsect}(end - nangovlap{radsect} + 1:end), size(tmp, 1), 1).*tmp(:, 1:nangovlap{radsect}, coil, slc) + ...
                        repmat(flipvec(ang_trans_wts{radsect}(1:nangovlap{radsect})), size(tmp1, 1), 1).*tmp1(:, end - nangovlap{radsect} + 1:end, coil, slc);
                    % Adjust end border.
                    tmp(:, end - nangovlap{radsect} + 1:end, coil, slc) = ...
                        repmat(flipvec(ang_trans_wts{radsect}(end - nangovlap{radsect} + 1:end)), size(tmp, 1), 1).*tmp(:, end - nangovlap{radsect} + 1:end, coil, slc) + ...
                        repmat(ang_trans_wts{radsect}(1:nangovlap{radsect}), size(tmp2, 1), 1).*tmp2(:, 1:nangovlap{radsect}, coil, slc);
                end
            end

            ang_range = [1:size(tmp, 2)] + (angsect - 1)*size(tmp, 2);
            xktmp{radsect}(:, ang_range, :, :) = tmp;
        end
    end
    %% Then combine and smooth borders of radial sectors.
    if nradsects == 1
        xk = xktmp{1};
    else
        xk = zeros(size(kdss_interp)); % [nx/2, ninterp, nvcoils, length(slctd)]
        for radsect = 1:nradsects
            if radsect == 1
                tmp = xktmp{radsect}(1:end - nradovlap{radsect}, :,:,:);

                tmp2 = xktmp{radsect + 1}; % Next sector.
                for coil = 1:nvcoils
                    for slc = 1:length(slctd)
                        % Adjust end border.
                        tmp(end - nradovlap{radsect} + 1:end, :, coil, slc) = ...
                            repmat(flipvec(rad_trans_wts{radsect}(end - nradovlap{radsect} + 1:end)).', 1, size(tmp, 2)).*tmp(end - nradovlap{radsect} + 1:end, :, coil, slc) + ...
                            repmat(rad_trans_wts{radsect}(1:nradovlap{radsect}).', 1, size(tmp2, 2)).*tmp2(1:nradovlap{radsect}, :, coil, slc);
                    end
                end
            elseif radsect == nradsects
                tmp = xktmp{radsect}(nradovlap{radsect} + 1:end, :,:,:);

                tmp1 = xktmp{radsect - 1}; % Previous sector.
                for coil = 1:nvcoils
                    for slc = 1:length(slctd)
                        % Adjust beginning border.
                        tmp(1:nradovlap{radsect}, :, coil, slc) = ...
                            repmat(rad_trans_wts{radsect}(end - nradovlap{radsect} + 1:end).', 1, size(tmp, 2)).*tmp(1:nradovlap{radsect}, :, coil, slc) + ...
                            repmat(flipvec(rad_trans_wts{radsect}(1:nradovlap{radsect})).', 1, size(tmp1, 2)).*tmp1(end - nradovlap{radsect} + 1:end, :, coil, slc);
                    end
                end
            else % Don't need to worry about boundaries.
                tmp = xktmp{radsect}(nradovlap{radsect} + 1:end - nradovlap{radsect}, :,:,:);

                tmp1 = xktmp{radsect - 1}; % Previous sector.
                tmp2 = xktmp{radsect + 1}; % Next sector.
                for coil = 1:nvcoils
                    for slc = 1:length(slctd)
                        % Adjust beginning border.
                        tmp(1:nradovlap{radsect}, :, coil, slc) = ...
                            repmat(rad_trans_wts{radsect}(end - nradovlap{radsect} + 1:end).', 1, size(tmp, 2)).*tmp(1:nradovlap{radsect}, :, coil, slc) + ...
                            repmat(flipvec(rad_trans_wts{radsect}(1:nradovlap{radsect})).', 1, size(tmp1, 2)).*tmp1(end - nradovlap{radsect} + 1:end, :, coil, slc);
                        % Adjust end border.
                        tmp(end - nradovlap{radsect} + 1:end, :, coil, slc) = ...
                            repmat(flipvec(rad_trans_wts{radsect}(end - nradovlap{radsect} + 1:end)).', 1, size(tmp, 2)).*tmp(end - nradovlap{radsect} + 1:end, :, coil, slc) + ...
                            repmat(rad_trans_wts{radsect}(1:nradovlap{radsect}).', 1, size(tmp2, 2)).*tmp2(1:nradovlap{radsect}, :, coil, slc);
                    end
                end
            end
            rad_range = [1:size(tmp, 1)] + sum(nrings(1:radsect - 1));
            xk(rad_range, :, :, :) = tmp;
        end
    end


    %% Transform single-slice k-space data into image data. Don't forget about the center k-space sample.
    x = zeros(nx, ny, nvcoils, length(slctd));
    if first_tp
        if do_test_trans
            xtest = zeros(nx, ny, nvcoils, length(slctd)); % Transform single-slice k-space data used for computation of weights, just as a test of the transform operation.
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
    end
    xk_demod = zeros(size(xk, 1)*size(xk, 2) + 1, nvcoils, length(slctd));
    if first_tp
        kdss_interp_demod = zeros(size(kdss_interp, 1)*size(kdss_interp, 2) + 1, nvcoils, length(slctd));
    end
    for slc = 1:length(slctd)
        for coil = 1:nvcoils
            xk_demod(:, coil, slc) = ...
                [col(xk(:,:, coil, slc).'); xk_cent(coil, slc)] .* ...
                exp(i * -1 * angle([col(mkss_interp(:,:, coil, slc).'); mkss(end, coil, slc)]));
            
            if first_tp
                kdss_interp_demod(:, coil, slc) = ...
                    [col(kdss_interp(:,:, coil, slc).'); kdss(end, coil, slc)] .* ...
                    exp(i * -1 * angle([col(mkss_interp(:,:, coil, slc).'); mkss(end, coil, slc)]));
            end
        end
    end
    dummy = rec_nufft(xk_demod, kinfo_interp, 'scaninfo', scaninfo, 'imsize', [nx ny], 'fm', fm, 'beta', 273, 'stop_diff_tol', 0, 'niter', 5);
    x = permute(dummy, [1 2 4 3]);
    if first_tp & do_test_trans
        dummy = rec_nufft(kdss_interp_demod, kinfo_interp, 'scaninfo', scaninfo, 'imsize', [nx ny], 'fm', fm, 'beta', 273, 'stop_diff_tol', 0, 'niter', 5);
        xtest = permute(dummy, [1 2 4 3]);
    end

    %% Combine image data from all coils into one final image.
    %xc(:,:,:, tp) = squeeze(sqrt(sum(abs(x).^2, 3))); % x combined using sum-of-squares.
    xc = squeeze(sqrt(sum(abs(x).^2, 3))); % x combined using sum-of-squares.
    %xc = coil_comb(x, c(:,:, slctd, :,:));
    %xc(:,:, tp) = x(:,:, 4); % Special case: Save just the 4th coil of 1 separated slice.
    %xc(:,:,:, tp) = x; % Special case: Save all coils of 1 separated slice.
    if first_tp & do_test_trans
        xtestc = squeeze(sqrt(sum(abs(xtest).^2, 3)));
    end

    %% Plots.
    if do_plots
        figure;
        %img(xc(:,:,:, tp), 'clims', [0 215]); % For CG with NUFFT. clim good for ACR phantom.
        img(xc(:,:,:, tp), 'clims', [0 65]); % For CG with NUFFT. clim good for human brain.
        %img(xc(:,:,:, tp), 'clims', [0 5.8e4]); % For gridding. clim good for ACR phantom.
        %img(xc(:,:,:, tp), 'clims', [0 2.1e4]); % For gridding. clim good for human brain.
        title(['GRAPPA SMS Recon: tp ' num2str(tp) ...
               ', ykern.coord.1: [' num2str(ykern_coord{1}) ']' ...
               ', ykern.coord.2: [' num2str(ykern_coord{2}) ']' ...
               ', ninterp: ' num2str(ninterp) ...
               ', nradsects: ' num2str(nradsects) ...
               ', nangsects: ' num2str(nangsects) ...
               ', nave: ' num2str(nave)])
        if first_tp & do_test_trans
            figure;
            %img(xtestc, 'clims', [0 215]); title('Testing transformation from k-space to image domain') % For CG with NUFFT. clim good for ACR phantom.
            img(xtestc, 'clims', [0 65]); title('Testing transformation from k-space to image domain') % For CG with NUFFT. clim good for human brain.
            %img(xtestc, 'clims', [0 5.8e4]); title('Testing transformation from k-space to image domain') % For gridding. clim good for ACR phantom.
            %img(xtestc, 'clims', [0 2.1e4]); title('Testing transformation from k-space to image domain') % For gridding. clim good for human brain.
        end
    end
end % if do_recon == 1

first_tp = 0;
