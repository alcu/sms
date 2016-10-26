function xk = grappa_recon_center(sig, w, nslcs)
% Uses GRAPPA weights for simultaneous multislice reconstruction. Produces
% k-space data for each coil and each separated slice for the center k-space
% sample only.
%
% Inputs:
% sig:   K-space data from simultaneous slices. [nsrc, ncoils, nacqs]
% w:     GRAPPA weights. [ncoils*nsrc, ntcoils*nslcs, nacqs]
% nslcs: number of simultaneous slices.
%
% Outputs:
% xk: Separated k-space data. [ntcoils, nslcs_overall]

%nslcs_overall = size(sigss, 4); % Overall number of slices (after separation).
nacqs = size(sig, 3); % Number of data acquisitions for this volume.
nslcs_overall = nslcs*nacqs;
%nslcs = nslcs_overall/nacqs; % Number of simultaneously acquired slices.
ncoils = size(sig, 2);
ntcoils = size(w, 2)/nslcs; % Number of coils in target.
nsrc = size(sig, 1);
%nx = size(sig, 1);
%ny = size(sig, 2);

xk = zeros(ntcoils, nslcs_overall);
for acq = 1:nacqs
    
    % Actual slices we are reconstructing.
    slcs_actual = [1 : nacqs : 1 + nacqs*(nslcs - 1)] + acq - 1;
    
    xk(:, slcs_actual) = reshape(reshape(sig(:, :, acq), 1, []) * squeeze(w(:, :, acq)), ntcoils, nslcs);
end
