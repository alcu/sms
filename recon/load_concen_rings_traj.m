%%% Concentric rings kx-ky trajectory.
load ../data/concen_ring/concen_rings2.mat % Vars: gx, gy, gz, gxy_beg_idx, kx_circ, ky_circ, kz_circ, kx_act, ky_act, kz_act, rast, fovk_unfudged, idx_circ. Note: These go from out to in. Careful about the units!
tmp_kx = kx_circ; % Need to swap kx and ky (and flip them) for image display to be correct.
kx_circ = -1*ky_circ(1:end - 1); % Delete 1 sample b/c scanner makes grad res divisible by 4 for some reason.
ky_circ = -1*tmp_kx(1:end - 1);
kz_circ = kz_circ(1:end - 1);
idx_circ = idx_circ(1:end - 1);
kinfo.k = nx/2*(kx_circ + i*ky_circ)/fovk_unfudged;
kinfo.t = [0 : 1 : idx_circ(end) - 1]*rast*1e-3;
kinfo.t = flipvec(scaninfo.te*1e-3 - kinfo.t); % Needs scaninfo from info.mat.
kinfo.t = kinfo.t(idx_circ);
scaninfo.revflg = 3; % Use tts = kinfo.t in rec_nufft(). Also needed for rec_fmap(), since it uses rec_nufft().
kinfo.kmod = ones(size(kinfo.k));
tmp_kx = kx_act;
kx_act = -1*ky_act(1:end - 1);
ky_act = -1*tmp_kx(1:end - 1);
kz_act = kz_act(1:end - 1);

% Density compensation factor.
k_act = nx/2*(kx_act + i*ky_act)/fovk_unfudged;
vd = 1; % Not variable density.
g = [diff(k_act) 0]; % Zero at end b/c out to in trajectory.
tmp_ang = unwrap(angle(kinfo.k(1:5)));
if tmp_ang(1) > tmp_ang(end) % Clockwise trajectory.
    kinfo.kdens = -1*abs(g).*sin(angle(g)-angle(k_act)).*vd;
else % Counterclockwise trajectory.
    kinfo.kdens = abs(g).*sin(angle(g)-angle(k_act)).*vd;
end
%kinfo.kdens = kinfo.kdens.*kinfo.kmod;
kinfo.kdens = kinfo.kdens(idx_circ);
clear vd g tmp_kx tmp_ang
