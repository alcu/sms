function w = grappa_weights_ave(sig, sigss, xkern_coord, ykern_coord, kern_cent, xtr, ytr, varargin)
% Computes GRAPPA weights for simultaneous multislice reconstruction. Uses
% multiple SMS volumes for computation to try to minimize effects from
% physiological variation.
%
% Inputs:
% sig:   Cell array of k-space data from simultaneous slices. Each cell: [krow, kcol, ncoils, nacqs]
% sigss: K-space data from separate slices (need to be Gz-modulated). [krow, kcol, ntcoils, nslcs_overall]
% xkern_coord: x kernel weight coordinates (can "skip over" certain source positions). [1, nrows]
% ykern_coord: y kernel weight coordinates (can "skip over" certain source positions). [1, ncols]
% kern_cent: Kernel center (coordinate for where the kernel computes the value). [1, 2]
% xtr: x training range. [1, num_x_training_pts] e.g. xtr = [kern_cent(1) : 2 : nx - kern_cent(1) + 1];
% ytr: y training range. [1, num_y_training_pts] e.g. ytr = [kern_cent(2) : 5 : ny - kern_cent(2) + 1];
% [trg_rad_offset]: radial offset for Strg (for asymmetric kernel weights).
%
% Outputs:
% w: GRAPPA weights. [ncoils*nsrc, ntcoils*nslcs, nacqs]

if nargin == 7
    trg_rad_offset = 0;
elseif nargin == 8
    trg_rad_offset = varargin{1};
else
    error('Number of args wrong.')
end

nslcs_overall = size(sigss, 4); % Overall number of slices (after separation).
nacqs = size(sig{1}, 4); % Number of data acquisitions for this volume.
nslcs = nslcs_overall/nacqs; % Number of simultaneously acquired slices.
ncoils = size(sig{1}, 3);
ntcoils = size(sigss, 3); % Number of coils in target.
nsrc = length(xkern_coord)*length(ykern_coord);
nx = size(sig{1}, 1);
ny = size(sig{1}, 2);
%xtr = [kern_cent(1):nx - kern_cent(1) + 1]; % x training range.
%ytr = [kern_cent(2):ny - kern_cent(2) + 1]; % y training range.
w = zeros(ncoils*nsrc, ntcoils*nslcs, nacqs);
nave = length(sig); % Number of SMS volumes to use for weight computation.

for acq = 1:nacqs
    
    % Actual slices we are reconstructing.
    slcs_actual = [1 : nacqs : 1 + nacqs*(nslcs - 1)] + acq - 1;
    
    Strg = sigss(xtr + trg_rad_offset, ytr, :, slcs_actual);
    Strg = reshape(Strg, [], ntcoils*nslcs);
    nrep = size(Strg, 1);
    Strg = repmat(Strg, nave, 1);
    Ssrc = zeros(nave*nrep, ncoils*nsrc);
    xpos = repmat([1:nx].', 1, ny); xpos = xpos(xtr, ytr); xpos = xpos(:);
    ypos = repmat([1:ny], nx, 1); ypos = ypos(xtr, ytr); ypos = ypos(:);
    for ii = 1:nave
        Ssrc_tmp = zeros(nrep, ncoils*nsrc);
        count = 1;
        for rep = 1:nrep
            xcur = xpos(count);
            ycur = ypos(count);
            % Compute coordinates of current kernel position.
            xkern_cur = xkern_coord + (xcur - kern_cent(1));
            ykern_cur = ykern_coord + (ycur - kern_cent(2));
            count = count + 1;
            for coil = 1:ncoils
                Ssrc_tmp(rep, (coil - 1)*nsrc + 1 : coil*nsrc) = reshape(sig{ii}(xkern_cur, ykern_cur, coil, acq), 1, []);
            end
        end
        Ssrc((ii - 1)*nrep + 1 : ii*nrep, :) = Ssrc_tmp;
    end
    w(:,:, acq) = pinv(Ssrc) * Strg;
    
end
