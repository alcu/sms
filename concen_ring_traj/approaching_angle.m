function out = approaching_angle(ref_angle, x, y)
% function out = approaching_angle(ref_angle, x, y)
%
% This determines whether the given [x y] vector is "approaching" the
% ref_angle, which is assumed to be a value from MATLAB's angle()
% function. i.e. If the [x y] vector is rotating clockwise, is it rotating
% toward or rotating away from the ref_angle. This is a bit tricky because of
% angle-wrapping of MATLAB's angle().

if ref_angle > 0
    if (angle(x + i*y) > 0) & (angle(x + i*y) > ref_angle)
        out = 1;
    elseif (angle(x + i*y) > 0) & (angle(x + i*y) < ref_angle)
        out = 0;
    elseif (angle(x + i*y) < 0) & (angle(x + i*y) < ref_angle - pi)
        out = 1;
    elseif (angle(x + i*y) < 0) & (angle(x + i*y) > ref_angle - pi)
        out = 0;
    elseif angle(x + i*y) == 0
        out = 0;
    elseif angle(x + i*y) == ref_angle - pi
        out = 0;
    elseif angle(x + i*y) == ref_angle
        out = 0;
    else
        error('Something is wrong...')
    end
elseif ref_angle < 0
    if (angle(x + i*y) < 0) & (angle(x + i*y) > ref_angle)
        out = 1;
    elseif (angle(x + i*y) < 0) & (angle(x + i*y) < ref_angle)
        out = 0;
    elseif (angle(x + i*y) > 0) & (angle(x + i*y) < ref_angle + pi)
        out = 1;
    elseif (angle(x + i*y) > 0) & (angle(x + i*y) > ref_angle + pi)
        out = 0;
    elseif angle(x + i*y) == 0
        out = 1;
    elseif angle(x + i*y) == ref_angle + pi
        out = 0;
    elseif angle(x + i*y) == ref_angle
        out = 0;
    else
        error('Something is wrong...')
    end
else % ref_angle == 0
    if angle(x + i*y) == pi
        out = 0;
    elseif angle(x + i*y) == 0
        out = 0;
    elseif angle(x + i*y) > 0
        out = 1;
    elseif angle(x + i*y) < 0
        out = 0;
    else
        error('Something is wrong...')
    end
end
