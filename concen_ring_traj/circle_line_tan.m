function [x y] = circle_line_tan(x1, y1, r1)
% function [x y] = circle_line_tan(x1, y1, r1)
% 
% Finds tangent point(s), of a circle with center (x1, y1) and radius r1, and
% a line with equation y = a*x. (Assumes possible tangent lines are NOT
% vertical. If any are, it throws an error.)
%
% Output:
%        - Point 1: x(1) = x-coord, y(1) = y-coord
%        - Point 2: x(2) = x-coord, y(2) = y-coord

b24ac = 64*x1^2*y1^2 - 4*(-4*x1^2 + 4*r1^2)*(-4*y1^2 + 4*r1^2);
if b24ac <= 0
    error('Circle is probably around the center point...')
elseif x1 == r1
    error('One possible tangent line is vertical: To do later...')
else
    a1 = (-8*x1*y1 + sqrt(b24ac))/(2*(-4*x1^2 + 4*r1^2));
    a2 = (-8*x1*y1 - sqrt(b24ac))/(2*(-4*x1^2 + 4*r1^2));
    
    % Point 1.
    x(1) = (2*x1 + 2*y1*a1)/(2*(1 + a1^2));
    y(1) = a1*x(1);
    
    % Point 2.
    x(2) = (2*x1 + 2*y1*a2)/(2*(1 + a2^2));
    y(2) = a2*x(2);
end

do_plots = 0;
if do_plots
    figure;
    theta = [0:0.01:2*pi];
    plot(x1 + r1*cos(theta), y1 + r1*sin(theta)); hold on
    plot([0 x(1)], [0 y(1)], '-r');
    plot([0 x(2)], [0 y(2)], '-m');
    plot(x(1), y(1), 'xg')
    plot(x(2), y(2), 'xk')
    legend('Circle', 'Line 1', 'Line 2', 'Tan pt 1', 'Tan pt 2')
    axis square
end
