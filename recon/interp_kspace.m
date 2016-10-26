function [kd_interp kl_interp tts_interp mk_interp] = interp_kspace(kd, kx_circ, ky_circ, tts, mk, ninterp, nacqs, ncoils, nx, varargin)
% Interpolates concentric ring k-space data to a constant angular velocity trajectory.
% 
% Inputs:
% kd:              K-space data.
% kx_circ:         Kx-space location.
% ky_circ:         Ky-space location.
% tts:             Vector of readout time. Enter as [] to not interpolate this.
% mk:              K-space modulation (from readout Gz). Enter as [] to not interpolate this.
% ninterp:         Number of evenly spaced samples to interpolate to.
% nacqs:           Number of acquisitions (or slices for single-slice data).
% ncoils:          Number of coils.
% nx:              Image matrix size.
% [ang_to_interp]: Angles to interpolate to. (IMPORTANT: Make sure angles are in decreasing order.)

if nargin == 9
    ang_to_interp = fliplr([-pi : 2*pi/ninterp : pi - 2*pi/ninterp]);
else
    ang_to_interp = varargin{1}(:).'; % Make into a row vector.
    ninterp = length(ang_to_interp);
end

% Find segment indices.
idx = zeros(nx/2, 2); idx(1, 1) = 1; % Index of starting, idx(seg_num, 1), and ending, idx(seg_num, 2), point of segment seg_num.
cur_rad = abs(kx_circ(idx(1, 1)) + i*ky_circ(idx(1, 1)));
seg_count = 1;
while true
    if cur_rad < 1.0 % Have reached the center of k-space, doesn't make sense to interpolate it.
        break;
    end
    
    idx_tmp = find(abs(abs(kx_circ(idx(seg_count, 1):end) + i*ky_circ(idx(seg_count, 1):end)) - cur_rad) > 1.0, ...
                   1, 'first');

    if isempty(idx_tmp) % This is the last segment.
        idx(seg_count, 2) = size(kd, 1);
        break;
    else
        idx(seg_count, 2) = idx(seg_count, 1) - 2 + idx_tmp;
    end

    seg_count = seg_count + 1;
    idx(seg_count, 1) = idx(seg_count - 1, 2) + 1;
    cur_rad = abs(kx_circ(idx(seg_count, 1)) + i*ky_circ(idx(seg_count, 1)));
end

kd_interp = zeros(nx/2, ninterp, ncoils, nacqs); % K-space data interpolated to constant ang velocity traj.
kl_interp = zeros(nx/2, ninterp); % K-space locations interpolated to constant ang velocity traj.
tts_interp = zeros(nx/2, ninterp); % Readout time vector interpolated to constant ang velocity traj.
mk_interp = zeros(nx/2, ninterp, ncoils, nacqs); % K-space modulation interpolated to constant ang velocity traj.
cur_klseg_angle = cell(1, nx/2); % Cell array of angle of k-space LOCATION segments, where each segment is along a ring (circle) with different radius.
cur_kseg_start_idx = zeros(1, nx/2);
ang_to_interp_shifted = zeros(nx/2, length(ang_to_interp));
for acq = 1:nacqs
    fprintf(['\nInterpolating k-space: acq ' num2str(acq) ' of ' num2str(nacqs) '...\n']);
    for coil = 1:ncoils
        for seg = 1:nx/2
            if (acq == 1) & (coil == 1)
                klseg = kx_circ(idx(seg, 1):idx(seg, 2)) + i*ky_circ(idx(seg, 1):idx(seg, 2));
                cur_klseg_angle{seg} = unwrap(angle(klseg));
                cur_kseg_start_idx(seg) = find(ang_to_interp < cur_klseg_angle{seg}(1), 1, 'first'); % Idx of first pt with angle less than cur_klseg_angle{seg}(1).
                ang_to_interp_shifted(seg, :) = unwrap(circshift(ang_to_interp, [0 -1*(cur_kseg_start_idx(seg) - 1)])); % ang_to_interp_shifted(seg, :)(1) is pt closest to (and <) cur_klseg_angle{seg}(1).
                kl_interp(seg, :) = circshift(interp1(cur_klseg_angle{seg}, klseg, ang_to_interp_shifted(seg, :), 'linear'), ...
                                              [0 (cur_kseg_start_idx(seg) - 1)]);

                if ~isempty(tts)
                    ttsseg = tts(idx(seg, 1):idx(seg, 2));
                    tts_interp(seg, :) = circshift(interp1(cur_klseg_angle{seg}, ttsseg, ang_to_interp_shifted(seg, :), 'linear'), ...
                                                   [0 (cur_kseg_start_idx(seg) - 1)]);
                end
            end

            kdseg = kd(idx(seg, 1):idx(seg, 2), coil, acq);
            kd_interp(seg, :, coil, acq) = circshift(interp1(cur_klseg_angle{seg}, kdseg, ang_to_interp_shifted(seg, :), 'linear'), ...
                                                     [0 (cur_kseg_start_idx(seg) - 1)]);

            if ~isempty(mk)
                mkseg = mk(idx(seg, 1):idx(seg, 2), coil, acq);
                mk_interp(seg, :, coil, acq) = circshift(interp1(cur_klseg_angle{seg}, mkseg, ang_to_interp_shifted(seg, :), 'linear'), ...
                                                         [0 (cur_kseg_start_idx(seg) - 1)]);
            end

            % Error checking for unwanted extrapolation.
            if sum(isnan(kd_interp(seg, :, coil, acq))) > 0
                error('Unwanted extrapolation has occured during interpolation.')
            end
        end
    end
end
