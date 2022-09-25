function [N,A] = trap(upper,lower,bigN,q_inf,p_inf,c)
% function for using composite trapezoidal rule
% inputs are the Cp distrubutions and # of panels to iterate over (bigN)
% outputs are the normal and axial force per unit span

yt = @(x) 12/100*(c/.2)*(.2969*sqrt(x/c)-.126*(x/c)-.3516*(x/c)^2+ ...
   .2843*(x/c)^3-.1036*(x/c)^4); % equation of thickness
start = 0; stop = c; % stop is = c = 2 m
X = linspace(start,stop,bigN+1);
N = 0; % preallocate Normal Force per Unit Span
A = 0; % preallocate Axial Force per Unit Span
for k=1:bigN % comp trap rule summation loop
    % Compute Upper and Lower Pressures at k and k+1 points 
    Cpl = fnval(lower,X(k)/c)*q_inf+p_inf; 
    Cpl2 = fnval(lower,X(k+1)/c)*q_inf+p_inf;
    Cpu = fnval(upper,X(k)/c)*q_inf+p_inf;
    Cpu2 = fnval(upper,X(k+1)/c)*q_inf+p_inf;
    % Compute changes in x and y
    dx = (X(k+1)-X(k));
    dy = (yt(X(k+1))-yt(X(k)));
    % Compute upper and lower suface contributions
    nu = -(Cpu2+Cpu)/2;
    nl = (Cpl2+Cpl)/2;
    % Compute upper and lower suface contributions
    au = (Cpu2+Cpu)/2;
    al = (Cpl2+Cpl)/2;
    % Sum
    N = N+dx*(nl+nu);
    A = A+dy*(al+au);% symmetrical airfoil
end

end