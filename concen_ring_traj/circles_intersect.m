function [x y] = circles_intersect(x1, y1, r1, x2, y2, r2, varargin)
% function [x y] = circles_intersect(x1, y1, r1, x2, y2, r2, [do_plots])
% 
% Finds intersection point(s), if any, of 2 circles, one with center (x1, y1)
% and radius r1, and the other with center (x2, y2) and center r2.
%
% Output:
%        If 2 intersection points, each of x and y will be of size [1 2]:
%        - Point 1: x(1) = x-coord, y(1) = y-coord
%        - Point 2: x(2) = x-coord, y(2) = y-coord
%        If 1 intersection point, each of x and y will be of size [1 1]:
%        - x = x-coord, y = y-coord
%        If no intersection points, x == y == [].

if x1 == x2 & y1 == y2 & r1 == r2
    error('You used the same circle for both circle inputs: Infinite number of intersections.')
    return
end

%do_plots = 1;
if nargin == 7
    do_plots = logical(varargin{1});
else
    do_plots = false;
end
if do_plots
    figure;
    theta = [0:0.01:2*pi];
    plot(x1 + r1*cos(theta), y1 + r1*sin(theta)); hold on
    plot(x2 + r2*cos(theta), y2 + r2*sin(theta), '-r');
    legend('Circle 1', 'Circle 2'); axis square
end

if y1 == y2
    if x1 == x2 % Circles have the same center.
        x = []; y = []; return
    else
        x(1) = ((x2^2 - x1^2) + (r1^2 - r2^2))/(2*(x2 - x1)); % Variable "c" in my notes.
        b24ac = 4*y1^2 - 4*(x(1)^2 - 2*x1*x(1) + x1^2 + y1^2 - r1^2);
        if b24ac < 0 % No intersection.
            x = []; y = []; return
        else
            y(1) = (2*y1 + sqrt(b24ac))/2;
            if b24ac == 0 % 1 point of intersection.
                return
            else % 2 points of intersection.
                x(2) = x(1);
                y(2) = (2*y1 - sqrt(b24ac))/2;
                return
            end
        end
    end
else
    a = (x2 - x1)/(y1 - y2); % Variable "a" in my notes.
    b = ((x1^2 - x2^2) + (y1^2 - y2^2) + (r2^2 - r1^2))/(2*(y1 - y2)); % Variable "b" in my notes.
    b24ac = (-2*x1 + 2*a*b - 2*y1*a)^2 - 4*(1 + a^2)*(x1^2 + b^2 - 2*y1*b + y1^2 - r1^2);
    if b24ac < 0 % No intersection.
        x = []; y = []; return
    else
        x(1) = ((2*x1 - 2*a*b + 2*y1*a) + sqrt(b24ac))/(2*(1 + a^2));
        y(1) = a*x(1) + b;
        if b24ac == 0 % 1 point of intersection.
            return
        else % 2 points of intersection.
            x(2) = ((2*x1 - 2*a*b + 2*y1*a) - sqrt(b24ac))/(2*(1 + a^2));
            y(2) = a*x(2) + b;
            return
        end
    end
end
