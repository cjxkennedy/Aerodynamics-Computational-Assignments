function [x,y,aLo,dCLda] = NACA_Airfoils(m,p,t,c,N)
%NACA_Airfoils func computes boundary points of airfoil, zero lift angle of
%attack and lift slope.
%
%   Uses vortex panel method to find characteristics of given NACA airfoil.
%   Inputs are parameters of NACA airfoils (maximum chamber, location of
%   max chamber, and thickness), chord length, and number of discrete
%   vortices, N.
%
%   Requires: Vortex_Panel.m to calculate coefficient of lift
%
%   Author: CJ Kennedy
%   Date: 4/3/21
%

thick = t; % redefine thickness :)
t = linspace(0,2*pi,N); % theta
x = c/2*(1+cos(t)); % transformation
yc = zeros(1,N); dycdx = zeros(1,N); % preallocate arrays

for i=1:N % loop thru N for yc and its deriv
    if x(i) < p*c % condition 1
         yc(i) = (m*x(i)/p^2)*(2*p-x(i)/c);
         dycdx(i) = (m*(2*p - x(i)/c))/p^2 - (m*x(i))/(c*p^2);
    else % condition 2
         yc(i) = m*(c-x(i))/(1-p)^2*(1+x(i)/c-2*p);
         dycdx(i) = (m*(c - x(i)))/(c*(p - 1)^2) - ...
            (m*(x(i)/c - 2*p + 1))/(p - 1)^2;
    end
end

Xi = atan(dycdx); % radians
yt = thick*(c/.2)*(.2969*sqrt(x./c)-.126*(x./c)-.3516*(x./c).^2+ ...
.2843*(x./c).^3-.1036*(x./c).^4); %thickness
x_R = zeros(1,N); y_R = zeros(1,N); % preallocate x and y 

for i=1:N % loop finding x and y location for N points
    if t(i) < pi % condition 1
        x_R(i) = x(i) + yt(i)*sin(Xi(i));
        y_R(i) = yc(i) - yt(i)*cos(Xi(i));
    else % condition 2
        x_R(i) = x(i) - yt(i)*sin(Xi(i));
        y_R(i) = yc(i) + yt(i)*cos(Xi(i));
    end
end

x = x_R; % reallocate x and y for function return
y = y_R; 

V_inf = 42.184; %Vinf from part 2 (m/s)
alpha = linspace(-15,15,N); % vary angle of attack (deg)
c_l = zeros(1,N); %preallocate c_l

for i=1:N % vortex loop to compute c_l
    [c_l(i), ~, ~, ~, ~] = Vortex_Panel(x,y,V_inf,alpha(i),0);
end
p = polyfit(alpha,c_l,1); % linear fit of aoa and cl
aLo = roots(p); % find x-int of fit
dCLda = (p(1)); % find first coeff / slope 
end