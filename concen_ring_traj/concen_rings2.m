%%% Generates k-space trajectory for concentric rings in kx-ky plane.
%%% Goes from inner ring to outer ring.
%%% Circles use maximum slew rate. Inspired by Jim Pipe's paper (MRM 2013).

% All units are in kHz, msec, mT, and m.
write = 0;                     % Write out gradient waveforms to disk.
write_dir = '~/write_dir/';    % Directory to write to. Don't forget trailing slash.
do_plots = 1;                  % Plot figures or not.
maxarray = 999999;             % Length of kx, ky, gx, and gy vectors.
matrix_size = 64;              % Images are matrix_size-by-matrix_size, and k-space will have matrix_size/2 circles (and a center point).
fovxy = 0.22;                  % (m)
resxy = fovxy/matrix_size;     % (m)
fovk_fudge = 1e-8;             % (cycles/m) Fudge factor to decrease fovk.
fovk = 1/resxy/2 - fovk_fudge; % Maximum extent of k-space. (cycles/m)
fovk_unfudged = 1/resxy/2;     % Maximum extent of k-space. (cycles/m)
slew_fudge = 0.01;             % (cycles/m) Fudge factor to decrease slew rate.
slewmax = 150 - slew_fudge;    % (mT/m/msec)
slewmax_unfudged = 150;        % (mT/m/msec)
gmax = 40;                     % (mT/m)
rast = 0.004;                  % Gradient raster time. (ms)
dgc = slewmax*rast;            % The most the gradients can change in 1 raster period. (mT/m)
gamma = 42.576;                % Typically 42.576 kHz/mT.
gamrast = gamma*rast;          % (cycles/mT)
first_gzblip_len = 18;         % Length (number of samples) of the first Gz blip (for SMS).
neg_gzblip_len = 26;           % Length (number of samples) of negative Gz blips (for SMS).
pos_gzblip_len = 18;           % Length (number of samples) of positive Gz blips (for SMS).
nslcs = 3;                     % Number of simultaneous slices for each acquisition.
dist_btw_slcs = 0.039;         % Distance between each of the simultaneous slices. (m)

% K-space trajectory. (cycles/m)
kx = zeros(1, maxarray);
ky = zeros(1, maxarray);
gx = zeros(1, maxarray); % (mT/m)
gy = zeros(1, maxarray); % (mT/m)

% Radii of concentric circles. (cycles/m)
rads = [fovk/(matrix_size/2) : fovk/(matrix_size/2) : fovk];

% Start out spiral going radially at max slew rate for 2 time points.
kx(1) = 0;
ky(1) = 0;
kx(2) = gamrast*dgc;
ky(2) = 0;
kx(3) = 3*gamrast*dgc;
ky(3) = 0;
gx(1) = 0;
gy(1) = 0;
gx(2) = (kx(2) - kx(1))/gamrast;
gy(2) = (ky(2) - ky(1))/gamrast;
gx(3) = (kx(3) - kx(2))/gamrast;
gy(3) = (ky(3) - ky(2))/gamrast;

idx = 4; % Index into kx, ky, gx, and gy.
idx_rad = 1; % Index into rads.
r = rads(idx_rad);
do_circle = 0;
slow_down = zeros(size(kx)); % Slow down trajectory for current sample.
do_first_trans = 1;
do_2nd_trans = 0;
found_centers = 0;
recorded_first_circ = 0;
recorded_first_trans1 = 0;
recorded_first_trans2 = 0;
do_break = 0;
trans_cutoffs = zeros(size(kx)); % Indicates where the circle transitions begin and end.
trans_cutoffs(2) = 1; % First transition (from center k-space point) begins at the 2nd sample.
ntrans = 0; % Number of transitions between circles.
done = 0;
while ~done
    %%% DEBUG
    % if idx == 63
    %     keyboard
    % end
    %%% END DEBUG
    ktx = kx(idx - 1) + gamrast*gx(idx - 1); % k point if keep same gradient.
    kty = ky(idx - 1) + gamrast*gy(idx - 1); % k point if keep same gradient.
    if do_circle % Have reached the next circle.
        %%% DEBUG
        % if idx_rad == 4
        %     keyboard
        % end
        %%% END DEBUG
        [kxtmp kytmp] = circles_intersect(ktx, kty, gamrast*dgc, 0, 0, r);
        if slow_down(idx)
            [kx(idx) ky(idx)] = pick_point_counterclock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
        else
            [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
        end
        gx(idx) = (kx(idx) - kx(idx - 1))/gamrast;
        gy(idx) = (ky(idx) - ky(idx - 1))/gamrast;
        if ~recorded_first_circ
            first_kxpoint = kx(idx); % First point on circle traj.
            first_kypoint = ky(idx);
            first_kpoint_angle = angle(kx(idx) + i*ky(idx));
            first_kpoint_idx = idx;
            trans_cutoffs(idx) = -1;
            recorded_first_circ = 1;
        end
        idx = idx + 1;
        do_checking = 0;
        while true
            %%% DEBUG
            % if idx == 140
            %     keyboard
            % end
            %%% END DEBUG
            if idx - 1 > first_kpoint_idx
                %%% DEBUG
                % if idx == 408 %idx_rad == 5 & angle(kx(idx - 1) + i*ky(idx - 1)) < 10*pi/180
                %     keyboard
                % end
                %%% END DEBUG
                if (~do_checking) & approaching_angle(first_kpoint_angle, kx(idx - 1), ky(idx - 1))
                    do_checking = 1;
                end
                if do_checking
                    [kxcompare kycompare] = pick_point_counterclock([kx(idx - 1) first_kxpoint], [ky(idx - 1) first_kypoint]);
                    if (kxcompare == first_kxpoint) & (kycompare == first_kypoint) % Completed circle trajectory. Go back to transitioning to the next circle.
                        kx(idx - 1) = 0; ky(idx - 1) = 0; % Erase previous point b/c it is beyond the first point.
                        if idx_rad == length(rads) % Done with entire trajectory.
                            done = 1;
                        else
                            idx = idx - 1; % Redo the erased point.
                            idx_rad = idx_rad + 1;
                            r = rads(idx_rad); % Radius of next circle.
                            do_circle = 0;
                        end
                        break
                    end
                end
            end
            ktx = kx(idx - 1) + gamrast*gx(idx - 1); % k point if keep same gradient.
            kty = ky(idx - 1) + gamrast*gy(idx - 1); % k point if keep same gradient.
            [kxtmp kytmp] = circles_intersect(ktx, kty, gamrast*dgc, 0, 0, r);
            if isempty(kxtmp) % K-space velocity is too fast to follow circle trajectory. Backup and slow down.
                %keyboard
                while slow_down(idx - 1) == 1
                    idx = idx - 1;
                end
                if idx <= first_kpoint_idx % Go back to the previous transition, and do it slower.
                    do_break = 1;
                    do_first_trans = 0;
                    do_2nd_trans = 1;
                    recorded_first_circ = 0;
                    trans_cutoffs(first_kpoint_idx) = 0;
                end
                if idx <= first_trans2_kpoint_idx % We have backed up far enough to get to the first transition.
                    do_first_trans = 1;
                    do_2nd_trans = 0;
                end
                if idx <= first_trans1_kpoint_idx % Backed up past the entire previous transition... Not sure what to do?
                    error('Backed up through an entire transition. Not sure what to do.')
                end
                slow_down(idx - 1) = 1;
                idx = idx - 1;
                if do_break % Go back to the previous transition, and do it slower.
                    do_break = 0;
                    do_circle = 0;
                    found_centers = 1; % Re-use previously found centers and radii for transitions.
                    break
                end
            else
                if slow_down(idx)
                    [kx(idx) ky(idx)] = pick_point_counterclock(kxtmp, kytmp);
                else
                    [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp);
                end
                gx(idx) = (kx(idx) - kx(idx - 1))/gamrast;
                gy(idx) = (ky(idx) - ky(idx - 1))/gamrast;
                idx = idx + 1;
            end
        end
    else % Do transitions, which are all semi-circles.
        if idx_rad == 1 % First transition from center to smallest circle.
            [kxtmp kytmp] = circles_intersect(ktx, kty, gamrast*dgc, 0, -r/2, r/2);
            if slow_down(idx)
                [kx(idx) ky(idx)] = pick_point_counterclock(kxtmp, kytmp);
                if abs(kx(idx) + i*ky(idx)) < abs(kx(idx - 1) + i*ky(idx - 1)) % If picking the counterclock pt results in going backward:
                    [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
                end
            else
                [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
            end
            gx(idx) = (kx(idx) - kx(idx - 1))/gamrast;
            gy(idx) = (ky(idx) - ky(idx - 1))/gamrast;
            if angle(kx(idx) + i*ky(idx)) <= -pi/2 % Done transitioning to first circle.
                if idx - 2 < first_gzblip_len % Trans is not long enough for first neg gz blip.
                    while slow_down(idx - 1) == 1
                        idx = idx - 1;
                    end
                    slow_down(idx - 1) = 1;
                    idx = idx - 1;
                else % Trans is long enough. Move on to the next circle.
                    do_circle = 1;
                    ntrans = ntrans + 1;
                    % Note: Don't change idx b/c want to overwrite this point during circle trajectory.
                end
            else
                idx = idx + 1;
            end
        else % Other transitions.
            if do_first_trans
                if ~recorded_first_trans1
                    first_trans1_kpoint_idx = idx; % Index of first point of first transition.
                    trans_cutoffs(idx) = 1;
                    recorded_first_trans1 = 1;
                    ntrans = ntrans + 1;
                end
                
                if ~found_centers
                    % Find center point of first transition circle.
                    if kx(idx - 1) == 0 % Vertical line.
                        kx_transcirc1_cent = 0;
                        ky_transcirc1_cent = r - 0.5*fovk/(matrix_size/2);
                    elseif ky(idx - 1) == 0 % Horizontal line.
                        kx_transcirc1_cent = r - 0.5*fovk/(matrix_size/2);
                        ky_transcirc1_cent = 0;
                    else
                        [kxtmp kytmp] = circle_line_intersect(0, 0, r - 0.5*fovk/(matrix_size/2), ky(idx - 1)/kx(idx - 1), 0, []);
                        if abs(kx(idx - 1) - kxtmp(1) + i*(ky(idx - 1) - kytmp(1))) < abs(kx(idx - 1) - kxtmp(2) + i*(ky(idx - 1) - kytmp(2))) % Pick point that's closer to [kx(idx - 1), ky(idx - 1)].
                            kx_transcirc1_cent = kxtmp(1);
                            ky_transcirc1_cent = kytmp(1);
                        else
                            kx_transcirc1_cent = kxtmp(2);
                            ky_transcirc1_cent = kytmp(2);
                        end
                    end
                    
                    % Find center point of 2nd transition circle.
                    [kxtmp kytmp] = circle_line_tan(kx_transcirc1_cent, ky_transcirc1_cent, 0.5*fovk/(matrix_size/2));
                    [kx_tan_pt ky_tan_pt] = pick_point_clock(kxtmp, kytmp);
                    ang_btw_cent_tan = abs(diff(unwrap([angle(kx_transcirc1_cent + i*ky_transcirc1_cent) angle(kx_tan_pt + i*ky_tan_pt)])));
                    polar_cent = abs(kx_transcirc1_cent + i*ky_transcirc1_cent)*exp(i*(angle(kx_tan_pt + i*ky_tan_pt) - ang_btw_cent_tan));
                    kx_transcirc2_cent = real(polar_cent);
                    ky_transcirc2_cent = imag(polar_cent);
                    
                    % Find radius of 2nd trans circle also.
                    rad_2nd_transcirc = 0.5*fovk/(matrix_size/2);
                    
                    found_centers = 1;
                end
                
                % Go along transition circle.
                [kxtmp kytmp] = circles_intersect(ktx, kty, gamrast*dgc, kx_transcirc1_cent, ky_transcirc1_cent, 0.5*fovk/(matrix_size/2));
                if isempty(kxtmp) % K-space velocity is too fast to follow circle trajectory. Backup and slow down.
                    %keyboard
                    while slow_down(idx - 1) == 1
                        idx = idx - 1;
                    end
                    if idx <= first_trans1_kpoint_idx % Go back to the previous circle, and do it slower.
                        do_circle = 1;
                        idx_rad = idx_rad - 1;
                        r = rads(idx_rad);
                        recorded_first_trans1 = 0;
                        recorded_first_circ = 1;
                        found_centers = 0;
                        trans_cutoffs(first_trans1_kpoint_idx) = 0;
                        ntrans = ntrans - 1;
                    end
                    slow_down(idx - 1) = 1;
                    idx = idx - 1;
                else
                    if slow_down(idx)
                        [kx(idx) ky(idx)] = pick_point_counterclock(kxtmp, kytmp);
                        if abs(kx(idx) + i*ky(idx)) < abs(kx(idx - 1) + i*ky(idx - 1)) % If picking the counterclock pt results in going backward:
                            [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
                        end
                    else
                        [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
                    end
                    gx(idx) = (kx(idx) - kx(idx - 1))/gamrast;
                    gy(idx) = (ky(idx) - ky(idx - 1))/gamrast;
                    if abs(kx(idx) + i*ky(idx)) > abs(kx_tan_pt + i*ky_tan_pt) % Start on next "transition circle."
                        idx = idx - 1; % Redo this point.
                        do_first_trans = 0;
                        do_2nd_trans = 1;
                    end
                    idx = idx + 1;
                end
            elseif do_2nd_trans % Do 2nd transition (using different "transition circle").
                %%% DEBUG
                % if idx == 77
                %     keyboard
                % end
                %%% END DEBUG
                if ~recorded_first_trans2
                    first_trans2_kpoint_idx = idx; % Index of first point of 2nd transition.
                    recorded_first_trans2 = 1;
                end
                
                % Go along transition circle.
                [kxtmp kytmp] = circles_intersect(ktx, kty, gamrast*dgc, kx_transcirc2_cent, ky_transcirc2_cent, rad_2nd_transcirc);
                if isempty(kxtmp) % K-space velocity is too fast to follow circle trajectory. Backup and slow down.
                    %keyboard
                    while slow_down(idx - 1) == 1
                        idx = idx - 1;
                    end
                    if idx <= first_trans2_kpoint_idx % Go back to the first transition, and do it slower.
                        do_first_trans = 1;
                        do_2nd_trans = 0;
                        recorded_first_trans2 = 0;
                    end
                    if idx <= first_trans1_kpoint_idx % We have backed up far enough to get to the previous circle.
                        do_circle = 1;
                        idx_rad = idx_rad - 1;
                        r = rads(idx_rad);
                        recorded_first_trans1 = 0;
                        recorded_first_circ = 1;
                        found_centers = 0;
                        trans_cutoffs(first_trans1_kpoint_idx) = 0;
                        ntrans = ntrans - 1;
                    end
                    slow_down(idx - 1) = 1;
                    idx = idx - 1;
                else
                    if slow_down(idx)
                        [kx(idx) ky(idx)] = pick_point_counterclock(kxtmp, kytmp);
                        if abs(kx(idx) + i*ky(idx)) < abs(kx(idx - 1) + i*ky(idx - 1)) % If picking the counterclock pt results in going backward:
                            [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
                        end
                    else
                        [kx(idx) ky(idx)] = pick_point_clock(kxtmp, kytmp); % Pick "most clockwise" point on circle.
                    end
                    gx(idx) = (kx(idx) - kx(idx - 1))/gamrast;
                    gy(idx) = (ky(idx) - ky(idx - 1))/gamrast;
                    if abs(abs(kx(idx) + i*ky(idx)) - r) < 0.03 % Done transitioning to circle.
                        %%% DEBUG
                        % if idx_rad == 3
                        %     keyboard
                        % end
                        %%% END DEBUG
                        if ((mod(ntrans, nslcs) == floor(nslcs/2) + 1) & (idx + 1 - first_trans1_kpoint_idx < neg_gzblip_len)) | ...
                           ((mod(ntrans, nslcs) ~= floor(nslcs/2) + 1) & (idx + 1 - first_trans1_kpoint_idx < pos_gzblip_len)) % Trans is not long enough for neg gz blip or for pos gz blip.
                            while slow_down(idx - 1) == 1
                                idx = idx - 1;
                            end
                            if idx <= first_trans2_kpoint_idx % Go back to the first transition, and do it slower.
                                do_first_trans = 1;
                                do_2nd_trans = 0;
                                recorded_first_trans2 = 0;
                            end
                            if idx <= first_trans1_kpoint_idx % We have backed up far enough to get to the previous circle.
                                do_circle = 1;
                                idx_rad = idx_rad - 1;
                                r = rads(idx_rad);
                                recorded_first_trans1 = 0;
                                recorded_first_circ = 1;
                                found_centers = 0;
                                trans_cutoffs(first_trans1_kpoint_idx) = 0;
                                ntrans = ntrans - 1;
                            end
                            slow_down(idx - 1) = 1;
                            idx = idx - 1;
                        else % Trans is long enough. Move on to the next circle.
                            do_circle = 1;
                            do_first_trans = 1;
                            do_2nd_trans = 0;
                            found_centers = 0;
                            recorded_first_trans1 = 0;
                            recorded_first_trans2 = 0;
                            recorded_first_circ = 0;
                        end
                    end
                    idx = idx + 1;
                end
            end
        end
    end
end


%%% Ramp gradients back to 0, then add appropriate gx and gy blips to get k-space traj back to 0.
zero_gx_idx = find(gx == 0, 5, 'first'); % Find indices of zero values in gx.
zero_gx_idx_diffs = diff(zero_gx_idx);
if zero_gx_idx_diffs(2:end) == ones(size(zero_gx_idx_diffs(2:end)))
    gxy_end_idx = zero_gx_idx(2) - 1; % Index of last point of "true" gx and gy waveforms.
    
    % Ramp gx to zero using max slew rate.
    ii = gxy_end_idx + 1;
    done_ramping_gx = 0;
    while ~done_ramping_gx
        if gx(gxy_end_idx) > 0
            gx(ii) = gx(ii - 1) - dgc;
        elseif gx(gxy_end_idx) < 0
            gx(ii) = gx(ii - 1) + dgc;
        end
        if ((gx(gxy_end_idx) > 0) & (gx(ii) < 0)) | ...
           ((gx(gxy_end_idx) < 0) & (gx(ii) > 0)) | ...
           (gx(ii) == 0)
            gx(ii) = 0;
            gx_len = ii; % Length of gx, including the ramp.
            done_ramping_gx = 1;
        else
            ii = ii + 1;
        end
    end
    
    % Ramp gy to zero using max slew rate.
    ii = gxy_end_idx + 1;
    done_ramping_gy = 0;
    while ~done_ramping_gy
        if gy(gxy_end_idx) > 0
            gy(ii) = gy(ii - 1) - dgc;
        elseif gy(gxy_end_idx) < 0
            gy(ii) = gy(ii - 1) + dgc;
        end
        if ((gy(gxy_end_idx) > 0) & (gy(ii) < 0)) | ...
           ((gy(gxy_end_idx) < 0) & (gy(ii) > 0)) | ...
           (gy(ii) == 0)
            gy(ii) = 0;
            gy_len = ii; % Length of gy, including the ramp.
            done_ramping_gy = 1;
        else
            ii = ii + 1;
        end
    end
else
    error('Have to manually check length of gx, gy waveforms.')
end
% Add gx and gy blips to get k-space traj back to 0.
gx_rewind_area = -1*sum(gx(1:gx_len))*rast; % (mT*ms/m)
gy_rewind_area = -1*sum(gy(1:gy_len))*rast; % (mT*ms/m)
gx_rewind_blip = 10*trapwave(gx_rewind_area/10/1000, rast*1e-3, gmax/10, slewmax/10*1000); % (mT/m)
gy_rewind_blip = 10*trapwave(gy_rewind_area/10/1000, rast*1e-3, gmax/10, slewmax/10*1000); % (mT/m)
gx = [gx(1:gx_len) gx_rewind_blip];
gy = [gy(1:gy_len) gy_rewind_blip];
% Make gx and gy same length.
if length(gx) > length(gy)
    gy = [gy zeros(1, length(gx) - length(gy))];
elseif length(gx) < length(gy)
    gx = [gx zeros(1, length(gy) - length(gx))];
end


% Compute lengths of transitions (in terms of number of samples).
trans_endpts = find(trans_cutoffs == -1);
trans_startpts = find(trans_cutoffs == 1);
trans_lengths = trans_endpts - trans_startpts;

% Plot figure of cutoff points and slow_down segments.
if do_plots
    figure;
    plot(trans_cutoffs)
    hold on;
    plot(slow_down, '--r')
    xlim([0 gxy_end_idx])
    legend('Transition cutoffs', 'Slow down segments')
end


%%% Create Gz waveform.
gz = zeros(size(gx));
dk = 1/(nslcs*dist_btw_slcs); % (cycles/m)
poswave = 10*trapwave(dk/(gamma*1000)/10, rast*1e-3, gmax/10, slewmax/10*1000); % (mT/m)
negwave = 10*trapwave(-1*(nslcs - 1)*dk/(gamma*1000)/10, rast*1e-3, gmax/10, slewmax/10*1000); % (mT/m)
if mod(matrix_size/2, nslcs) < floor(nslcs/2) + 1
    fact = -1*mod(matrix_size/2, nslcs);
elseif mod(matrix_size/2, nslcs) > floor(nslcs/2) + 1
    fact = nslcs - mod(matrix_size/2, nslcs);
elseif mod(matrix_size/2, nslcs) == floor(nslcs/2) + 1
    fact = round(nslcs/2) - 1;
end
lastwave = 10*trapwave(fact*dk/(gamma*1000)/10, rast*1e-3, gmax/10, slewmax/10*1000); % (mT/m)
blip_count = 1;
gz_cutoffs = zeros(size(gz)); % Cutoffs for gz blips.
for ii = 1:length(trans_cutoffs)
    if trans_cutoffs(ii) == 1
        gz_cutoffs(ii + 1) = 1;
        if mod(blip_count, nslcs) == floor(nslcs/2) + 1 % Neg blip.
            gz(ii + 1 : ii + length(negwave)) = negwave;
            gz_cutoffs(ii + length(negwave)) = -1;
        else % Pos blip.
            gz(ii + 1 : ii + length(poswave)) = poswave;
            gz_cutoffs(ii + length(poswave)) = -1;
        end
        blip_count = blip_count + 1;
    end
end
gz(end - length(lastwave) + 1 : end) = lastwave; % Last blip to rewind gz.
gz_cutoffs(end - length(lastwave) + 1) = 1;
gz_cutoffs(end) = -1;


%% Adjust transition cutoffs.
% Shift trans cutoffs over by 1 (b/c g is 1-off relative to k...)
ii = 1;
while ii < length(trans_cutoffs)
    if trans_cutoffs(ii) == 1
        trans_cutoffs(ii) = 0;
        trans_cutoffs(ii + 1) = 1;
        ii = ii + 2;
    elseif trans_cutoffs(ii) == -1
        % Do nothing, since the right edge of trans_cutoffs is actually the first point beyond the transition.
        ii = ii + 1;
    else
        ii = ii + 1;
    end
end
% Make same length as gx, gy, gz.
if length(trans_cutoffs) < length(gx)
    trans_cutoffs = [trans_cutoffs zeros(1, length(gx) - length(trans_cutoffs))];
elseif length(trans_cutoffs) > length(gx)
    trans_cutoffs = trans_cutoffs(1:length(gx));
end


%%% Reverse gx, gy, gz (to get a "circle-in").
gx = flipvec(gx);
gy = flipvec(gy);
gz = flipvec(gz);
gz_cutoffs = flipvec(gz_cutoffs);
trans_cutoffs = flipvec(trans_cutoffs);
%% Note: If playing gx and gy backwards (to get a "circle-in"), first point of
%% actual trajectory is: length(gx) - gxy_end_idx + 1.
%% Should begin the data readout at this point.
gxy_beg_idx = length(gx) - gxy_end_idx + 1;


%%% Compute actual kx, ky, kz values from gx, gy, gz.
gxc = cumsum(gx);
gyc = cumsum(gy);
gzc = cumsum(gz);
kx_act = gamrast*gxc; % Entire k-space trajectory (including preblips and transitions).
ky_act = gamrast*gyc;
kz_act = gamrast*gzc;
beg_trans_cut_idx = find(trans_cutoffs == -1); % Note switch (beg is end, and end is beg now).
end_trans_cut_idx = find(trans_cutoffs == 1);
% kx_circ = gamrast*gxc(gxy_beg_idx:beg_trans_cut_idx(1) - 1); % kx values for just the concentric circles (no transitions).
% ky_circ = gamrast*gyc(gxy_beg_idx:beg_trans_cut_idx(1) - 1); % ky values for just the concentric circles (no transitions).
% kz_circ = gamrast*gzc(gxy_beg_idx:beg_trans_cut_idx(1) - 1); % kz values for just the concentric circles (no transitions).
% idx_circ = [gxy_beg_idx:beg_trans_cut_idx(1) - 1]; % Indices of g[x,y,z]c for just the concentric circles (no transitions).
kx_circ = gamrast*gxc(gxy_beg_idx - 1:beg_trans_cut_idx(1)); % kx values for just the concentric circles (but with 1 transition pt on either side of the circle segment for interpolation during GRAPPA recon).
ky_circ = gamrast*gyc(gxy_beg_idx - 1:beg_trans_cut_idx(1)); % ky values for just the concentric circles (...).
kz_circ = gamrast*gzc(gxy_beg_idx - 1:beg_trans_cut_idx(1)); % kz values for just the concentric circles (...).
idx_circ = [gxy_beg_idx - 1:beg_trans_cut_idx(1)]; % Indices of g[x,y,z]c for just the concentric circles (...).
for ii = 2:length(beg_trans_cut_idx)
    % kx_circ = [kx_circ gamrast*gxc(end_trans_cut_idx(ii - 1) + 1:beg_trans_cut_idx(ii) - 1)];
    % ky_circ = [ky_circ gamrast*gyc(end_trans_cut_idx(ii - 1) + 1:beg_trans_cut_idx(ii) - 1)];
    % kz_circ = [kz_circ gamrast*gzc(end_trans_cut_idx(ii - 1) + 1:beg_trans_cut_idx(ii) - 1)];
    % idx_circ = [idx_circ [end_trans_cut_idx(ii - 1) + 1:beg_trans_cut_idx(ii) - 1]];
    kx_circ = [kx_circ gamrast*gxc(end_trans_cut_idx(ii - 1):beg_trans_cut_idx(ii))];
    ky_circ = [ky_circ gamrast*gyc(end_trans_cut_idx(ii - 1):beg_trans_cut_idx(ii))];
    kz_circ = [kz_circ gamrast*gzc(end_trans_cut_idx(ii - 1):beg_trans_cut_idx(ii))];
    idx_circ = [idx_circ [end_trans_cut_idx(ii - 1):beg_trans_cut_idx(ii)]];
end
% kx_circ = [kx_circ gamrast*gxc(end_trans_cut_idx(end) + 1:end)];
% ky_circ = [ky_circ gamrast*gyc(end_trans_cut_idx(end) + 1:end)];
% kz_circ = [kz_circ gamrast*gzc(end_trans_cut_idx(end) + 1:end)];
% idx_circ = [idx_circ [end_trans_cut_idx(end) + 1:length(gxc)]];
kx_circ = [kx_circ gamrast*gxc(end_trans_cut_idx(end):end)];
ky_circ = [ky_circ gamrast*gyc(end_trans_cut_idx(end):end)];
kz_circ = [kz_circ gamrast*gzc(end_trans_cut_idx(end):end)];
idx_circ = [idx_circ [end_trans_cut_idx(end):length(gxc)]];
if ~isequal(kx_circ, gamrast*gxc(idx_circ)) % Test if idx_circ is correct.
    error('idx_circ not correct')
end


%%% Write external waveforms.
if write
    % Write gradient waveforms to binary files. (*even* short integers -- the psd sets the EOS bit, so don't have to worry about it here.)
    maxiamp = 2^15 - 2;          % max instruction amplitude (max value of signed short)

    a_gx = max(abs(gx)); % (mT/m)
    tmp = 2*round((gx/a_gx)*maxiamp/2);
    tmp(end) = tmp(end) - 1;
    fid = fopen([write_dir 'gx.raw'], 'w', 'ieee-be');
    fwrite(fid, tmp, 'int16');
    fclose(fid);

    a_gy = max(abs(gy)); % (mT/m)
    tmp = 2*round((gy/a_gy)*maxiamp/2);
    tmp(end) = tmp(end) - 1;
    fid = fopen([write_dir 'gy.raw'], 'w', 'ieee-be');
    fwrite(fid, tmp, 'int16');
    fclose(fid);

    a_gzmod = max(abs(gz)); % (mT/m)
    tmp = 2*round((gz/a_gzmod)*maxiamp/2);
    tmp(end) = tmp(end) - 1;
    fid = fopen([write_dir 'gzmod.raw'], 'w', 'ieee-be');
    fwrite(fid, tmp, 'int16');
    fclose(fid);

    gres = length(gx);
    if gres ~= length(gy) | gres ~= length(gz)
        error('Something is wrong with length (res) of gradient waveforms.')
    end
    pw_gx = gres*rast; % (ms)
    
    %%
    %% Write info to text file.
    %%
    fid = fopen([write_dir 'multslc_arbtraj_info.txt'], 'w', 'ieee-be');
    fprintf(fid, 'float gxinfo[3] = {%.9f, %.9f, %.9f};\n\n', ...
            [a_gx/10; gres; pw_gx*1e3]);
    fprintf(fid, 'float gyinfo[1] = {%.9f};\n\n', ...
            a_gy/10);
    fprintf(fid, 'float gzinfo[1] = {%.9f};\n\n', ...
            a_gzmod/10);
    fclose(fid);
    
    
    % Reminder to run 'xlatebin' on the raw binary files.
    fprintf('*** Remember to run:\n"xlatebin -o gx.wav gx.raw; xlatebin -o gy.wav gy.raw; xlatebin -o gzmod.wav gzmod.raw"\n...to convert them to the format used by EXTWAVE().\n\n');
    
    % Reminder to update the .e file.
    fprintf('*** Also remember to update the gxinfo[], gyinfo[], and gzinfo[] arrays in the .e file.\n\n');
end


% Plot kx-ky traj from gx and gy.
if do_plots
    figure;
    plot(gamrast*cumsum(gx), gamrast*cumsum(gy), '-x'); axis square
    hold on
    for ii = 1:length(gz_cutoffs)
        if gz_cutoffs(ii) == 1
            plot(gamrast*gxc(ii), gamrast*gyc(ii), 'xr', 'linewidth', 2)
        elseif gz_cutoffs(ii) == -1
            plot(gamrast*gxc(ii), gamrast*gyc(ii), 'xg', 'linewidth', 2)
        end
    end
    for ii = 1:length(trans_cutoffs)
        if trans_cutoffs(ii) == 1
            plot(gamrast*gxc(ii), gamrast*gyc(ii), 'x', 'linewidth', 2)
        elseif trans_cutoffs(ii) == -1
            plot(gamrast*gxc(ii), gamrast*gyc(ii), 'xc', 'linewidth', 2)
        end
    end
    xlabel('k_x (cycles/m)'), ylabel('k_y (cycles/m)')

    figure;
    plot(kx_circ, ky_circ, '-x'); axis square
    xlabel('k_x (cycles/m)'), ylabel('k_y (cycles/m)')
    title('Non-transition samples only')

    % Plot 3D kx-ky-kz traj from gx, gy, and gz.
    figure;
    plot3(gamrast*cumsum(gx), gamrast*cumsum(gy), gamrast*cumsum(gz), '-'); axis square
    xlabel('k_x (cycles/m)'), ylabel('k_y (cycles/m)'), zlabel('k_z (cycles/m)')

    figure;
    plot3(kx_circ, ky_circ, kz_circ, '-'); axis square
    xlabel('k_x (cycles/m)'), ylabel('k_y (cycles/m)'), zlabel('k_z (cycles/m)')
    title('Non-transition samples only')

    % Plot gx, gy, gz.
    figure;
    plot(gx); hold on
    plot(gy, '-r')
    plot(gz, '-m')
    ylabel('(mT/m)'); xlabel('Sample number')
    legend('G_x', 'G_y', 'G_z')

    % Plot diff of gx, gy, gz.
    figure;
    plot(diff(gx)); hold on
    plot(diff(gy), '-r')
    plot(diff(gz), '-m')
    ylabel('(mT/m)')
    legend('Diff of G_x', 'Diff of G_y', 'Diff of G_z')

    % Plot kz.
    figure;
    plot(kz_act)
    ylabel('k_z (cycles/m)'), xlabel('Sample number')
    title('kz.act')
end
