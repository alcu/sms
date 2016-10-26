function [kdc A S] = coil_compr(kd, nvcoil, tp, varargin)
% function [kdc A S] = coil_compr(kd, nvcoil, tp, [A])
%
% This function does a simple SVD-based coil compression.
%
% Inputs:
% kd: K-space data to be compressed. [nksamp, ncoil, nacq, ntp]
% nvcoil: Desired number of virtual coils. [1, 1]
% tp: Vector of time points to compress. Use 0 to not compress anything. [1, ntp (or less)]
% [A]: Compression matrices to use. If not specified, the function computes
%      them and outputs them as A. [ncoil, nvcoil, nacq]
%
% Outputs:
% kdc: Compressed k-space data. [nksamp, nvcoil, nacq, length(tp)]
% A: Compression matrices. [ncoil, nvcoil, nacq]
% S: Vector of computed singular values if computing compression
%    matrices. [ncoil, nacq]

nksamp = size(kd, 1);
ncoil = size(kd, 2);
nacq = size(kd, 3);
ntp = size(kd, 4); % Number of time frames to use for coil compression matrix computation.

if nargin == 3
    comp_A = 1;
    A = zeros(ncoil, nvcoil, nacq);
    S = zeros(ncoil, nacq);
else
    comp_A = 0;
    A = varargin{1};
    if size(A, 3) ~= nacq
        error('Compression matrix input has a different number of acqs than data input.')
    end
    S = [];
end

if tp
    kdc = zeros(nksamp, nvcoil, nacq, length(tp));
else
    kdc = [];
end
for acq = 1:nacq
    if comp_A == 1
        % kdtmp = zeros(ntp*nksamp, ncoil);
        % for tp_idx = 1:ntp
        %     kdtmp([1:nksamp] + (tp_idx - 1)*nksamp, :) = kd(:,:, acq, tp_idx);
        % end
        kdtmp = reshape(permute(squeeze(kd(:,:, acq, :)), [1 3 2]), ntp*nksamp, ncoil);
        [U Stmp V] = svd(kdtmp, 'econ');
        S(:, acq) = diag(Stmp);
        A(:,:, acq) = V(:, 1:nvcoil);
    end
    
    if tp
        for tp_idx = 1:length(tp)
            tptd = tp(tp_idx);
            
            kdc(:,:, acq, tp_idx) = kd(:,:, acq, tptd)*A(:,:, acq);
        end
    end
end
