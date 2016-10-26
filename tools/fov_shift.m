function out = fov_shift(kd, kshmod)

ncoils = size(kd, 2);
nslcs = size(kd, 3);
out = zeros(size(kd));
for slc = 1:nslcs
    for coil = 1:ncoils
        out(:, coil, slc) = kshmod.*kd(:, coil, slc);
    end
end
