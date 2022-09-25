function [aLo] = aLoForThinCamber(m,p,c)
%aLoForThinCamber find the zero lift angle of attack for chambered airfoil
%using thin airfoil theory.
%
%   Computes proper integral over surfaces of airfoil. Inputs are NACA
%   airfoil parameters m and p, and chord length, c. Outputted zero-lift
%   angle of attack is in deg
%
%
%   Author: CJ Kennedy
%   Date: 4/3/21
%
syms t % symbolics theta, to be subsituted later
x = c/2*(1+cos(t)); % transformation
% Setup dzdx from NACA_Airfoils func
dycdxU = (m*(2*p - x/c))/p^2 - (m*x)/(c*p^2); 
dycdxL = (m*(c - x))/(c*(p - 1)^2) - (m*(x/c - 2*p + 1))/(p - 1)^2;
% Setup integrals to be evaluated
S1 = 1/pi*dycdxU*(cos(t)-1);
S2 = 1/pi*dycdxL*(cos(t)-1);
int1 = int(S1,t,0,1.77); % compute integral over first range
int2 = int(S2,t,1.77,pi); % compute integral over second range
aLo = rad2deg(double((int1+int2))); % return and convert to deg
end

