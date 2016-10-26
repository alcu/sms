function w = grappa_weights_center_ave(sig, sigss)
% Computes GRAPPA weights for simultaneous multislice reconstruction of the
% center k-space sample only. Uses multiple SMS volumes for computation to
% try to minimize effects from physiological variation.
%
% Inputs:
% sig:   Cell array of k-space data from simultaneous slices. Each cell: [nsrc, ncoils, nacqs]
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
nacqs = size(sig{1}, 3); % Number of data acquisitions for this volume.
nslcs = nslcs_overall/nacqs; % Number of simultaneously acquired slices.
ncoils = size(sig{1}, 2);
nsrc = size(sig{1}, 1);
%nx = size(sig, 1);
%ny = size(sig, 2);
%xtr = [2:nx - 1]; % x training range.
%ytr = [2:ny - 1]; % y training range.
w = zeros(ncoils*nsrc, ntcoils*nslcs, nacqs);
nave = length(sig); % Number of SMS volumes to use for weight computation.

for acq = 1:nacqs
    
    % Actual slices we are reconstructing.
    slcs_actual = [1 : nacqs : 1 + nacqs*(nslcs - 1)] + acq - 1;
    
    if isvector(sigss)
        Strg = sigss(slcs_actual);
    else
        Strg = sigss(:, slcs_actual);
    end
    Strg = reshape(Strg, [], ntcoils*nslcs);
    nrep = size(Strg, 1); % Should be 1.
    Strg = repmat(Strg, nave, 1);
    Ssrc = zeros(nave*nrep, ncoils*nsrc);
    %xpos = repmat([1:nx].', 1, ny); xpos = xpos(xtr, ytr); xpos = xpos(:);
    %ypos = repmat([1:ny], nx, 1); ypos = ypos(xtr, ytr); ypos = ypos(:);
    for ii = 1:nave
        Ssrc_tmp = zeros(nrep, ncoils*nsrc);
        %count = 1;
        for rep = 1:nrep
            %xcur = xpos(count);
            %ycur = ypos(count);
            %count = count + 1;
            for coil = 1:ncoils
                Ssrc_tmp(rep, (coil - 1)*nsrc + 1 : coil*nsrc) = reshape(sig{ii}(:, coil, acq), 1, []);
            end
        end
        Ssrc((ii - 1)*nrep + 1 : ii*nrep, :) = Ssrc_tmp;
    end
    w(:,:, acq) = pinv(Ssrc) * Strg;
    
end
