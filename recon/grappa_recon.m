function xk = grappa_recon(sig, w, nslcs, xkern_coord, ykern_coord, kern_cent, xtr, ytr)
% Uses GRAPPA weights for simultaneous multislice reconstruction. Produces
% k-space data for each coil and each separated slice.
%
% Inputs:
% sig:   K-space data from simultaneous slices. [krow, kcol, ncoils, nacqs]
% w:     GRAPPA weights. [ncoils*nsrc, ntcoils*nslcs, nacqs]
% nslcs: number of simultaneous slices.
% xkern_coord: x kernel weight coordinates (can "skip over" certain source positions). [1, nrows]
% ykern_coord: y kernel weight coordinates (can "skip over" certain source positions). [1, ncols]
% kern_cent: Kernel center (coordinate for where the kernel computes the value). [1, 2]
% xtr: x training range. [1, num_x_training_pts] e.g. xtr = [kern_cent(1) : 2 : nx - kern_cent(1) + 1];
% ytr: y training range. [1, num_y_training_pts] e.g. ytr = [kern_cent(2) : 5 : ny - kern_cent(2) + 1];
%
% Outputs:
% xk: Separated k-space data. [length(xtr), length(ytr), ntcoils, nslcs_overall]

%nslcs_overall = size(sigss, 4); % Overall number of slices (after separation).
nacqs = size(sig, 4); % Number of data acquisitions for this volume.
nslcs_overall = nslcs*nacqs;
%nslcs = nslcs_overall/nacqs; % Number of simultaneously acquired slices.
ncoils = size(sig, 3);
ntcoils = size(w, 2)/nslcs; % Number of coils in target.
nsrc = length(xkern_coord)*length(ykern_coord);
nx = size(sig, 1);
ny = size(sig, 2);
nrep = length(xtr)*length(ytr);

xk = zeros(length(xtr), length(ytr), ntcoils, nslcs_overall);
for acq = 1:nacqs
    
    % Actual slices we are reconstructing.
    slcs_actual = [1 : nacqs : 1 + nacqs*(nslcs - 1)] + acq - 1;
    
    % Using imfilter().
    %{
    for coil = 1:ncoils
        for slc = 1:nslcs
            for srccoil = 1:ncoils
                xk(:, :, coil, slcs_actual(slc)) = ...
                    xk(:, :, coil, slcs_actual(slc)) + ...
                    imfilter(sig(:,:, srccoil, acq), ...
                             reshape(w((srccoil - 1)*nsrc + 1:srccoil*nsrc, (slc - 1)*ncoils + coil, acq), ...
                                     sqrt(nsrc), sqrt(nsrc)), ...
                             'replicate');
            end
        end
    end
    %}
    
    % Using matrix multiplication.
    xpos = repmat([1:nx].', 1, ny); xpos = xpos(xtr, ytr); xpos = xpos(:);
    ypos = repmat([1:ny], nx, 1); ypos = ypos(xtr, ytr); ypos = ypos(:);
    Ssrc_tmp = zeros(nrep, ncoils*nsrc);
    for rep = 1:nrep
        xcur = xpos(rep);
        ycur = ypos(rep);
        % Compute coordinates of current kernel position.
        xkern_cur = xkern_coord + (xcur - kern_cent(1));
        ykern_cur = ykern_coord + (ycur - kern_cent(2));
        for coil = 1:ncoils
            Ssrc_tmp(rep, (coil - 1)*nsrc + 1 : coil*nsrc) = reshape(sig(xkern_cur, ykern_cur, coil, acq), 1, []);
        end
    end
    xk_tmp = Ssrc_tmp * w(:,:, acq);
    xk(:,:,:, slcs_actual) = reshape(xk_tmp, length(xtr), length(ytr), ntcoils, nslcs);
end
