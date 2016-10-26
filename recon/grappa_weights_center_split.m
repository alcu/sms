function w = grappa_weights_center_split(nslcs, sig, sigss)
% Computes GRAPPA weights for simultaneous multislice reconstruction of the
% center k-space sample only.
% Uses split slice-GRAPPA.
%
% Inputs:
% nslcs: Number of simultaneous slices. [1, 1]
% sig:   K-space data from simultaneous slices. [nsrc, ncoils, nacqs]
% sigss: Center k-space data from separate slices (need to be Gz-modulated). [ntcoils, nslcs_overall]
%
% Outputs:
% w:     GRAPPA weights. [ncoils*nsrc, ntcoils*nslcs, nacqs]

if isvector(sigss)
    nslcs_overall = length(sigss); % Special case for 1 target coil.
    ntcoils = 1; % Number of coils in target.
else
    nslcs_overall = size(sigss, 2); % Overall number of slices (after separation).
    ntcoils = size(sigss, 1); % Number of coils in target.
end
nacqs = nslcs_overall/nslcs; % Number of data acquisitions for this volume.
ncoils = size(sig, 2);
nsrc = size(sig, 1);
%nx = size(sig, 1);
%ny = size(sig, 2);
%xtr = [2:nx - 1]; % x training range.
%ytr = [2:ny - 1]; % y training range.
w = zeros(ncoils*nsrc, ntcoils*nslcs, nacqs);

for acq = 1:nacqs
    
    % Actual slices we are reconstructing.
    slcs_actual = [1 : nacqs : 1 + nacqs*(nslcs - 1)] + acq - 1;
    
    if isvector(sigss)
        Strg_tmp = sigss(slcs_actual);
    else
        Strg_tmp = sigss(:, slcs_actual);
    end
    nrep = 1; % Because just 1 target (the center).
    Strg = zeros(nrep*nslcs, nslcs*ntcoils);
    if isvector(sigss)
        for slc = 1:nslcs
            Strg([1:nrep] + (slc - 1)*nrep, [1:ntcoils] + (slc - 1)*ntcoils) = ...
                reshape(Strg_tmp(slc), nrep, ntcoils);
        end
    else
        for slc = 1:nslcs
            Strg([1:nrep] + (slc - 1)*nrep, [1:ntcoils] + (slc - 1)*ntcoils) = ...
                reshape(Strg_tmp(:, slc), nrep, ntcoils);
        end
    end
    Ssrc = zeros(nslcs*nrep, ncoils*nsrc);
    %xpos = repmat([1:nx].', 1, ny); xpos = xpos(xtr, ytr); xpos = xpos(:);
    %ypos = repmat([1:ny], nx, 1); ypos = ypos(xtr, ytr); ypos = ypos(:);
    sig_cur = sig(:,:, slcs_actual);
    for slc = 1:nslcs
        Ssrc_tmp = zeros(nrep, ncoils*nsrc);
        %count = 1;
        for rep = 1:nrep
            %xcur = xpos(count);
            %ycur = ypos(count);
            %count = count + 1;
            for coil = 1:ncoils
                Ssrc_tmp(rep, (coil - 1)*nsrc + 1 : coil*nsrc) = reshape(sig_cur(:, coil, slc), 1, []);
            end
        end
        Ssrc((slc - 1)*nrep + 1 : slc*nrep, :) = Ssrc_tmp;
    end
    w(:,:, acq) = pinv(Ssrc) * Strg;
    
end
