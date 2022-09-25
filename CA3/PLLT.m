function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%PLLT func computes span eff. factor, coefficient of Lift and Induced Drag
%
%   Uses fundamental eqn. of Prandtl Lifting Line Theory for finite wings
%   with thick airfoils. Allows for linear spanwise variation of inputs. 
%   Inputs are physical parameters of airfoil and N, number of odd terms in
%   expansion.
%
%
%   Author: CJ Kennedy
%   Date: 4/3/21
%

S = b/2*(c_r+c_t); % span calculations
AR = b^2/S; % aspect ratio

t = pi/(2*N):pi/(2*N):pi/2; % theta span
y = b/2*cos(t); % transformation

% define variables w linear spanwise variation across chord
c = c_r+(c_t-c_r)/(b/2)*y; 
a0 = a0_r+(a0_t-a0_r)/(b/2)*y;
geo = geo_r+(geo_t-geo_r)/(b/2)*y;
aero = aero_r+(aero_t-aero_r)/(b/2)*y;

% setup linear matrices in form Ax = b
smallB = zeros(1,N); bigA = zeros(N,N);
for i=1:N % go through columns
    ind = -1; % reset index
    for j=1:N % go through rows
        ind = ind+2; % increment (2*j-1)
        bigA(i,j) = 4*b/(a0(i)*c(i))*sin(ind*t(i))+ ...
            ind*sin(ind*t(i))/sin(t(i)); % A
    end 
    smallB(i) = geo(i)-aero(i); % b
end

x = linsolve(bigA,smallB'); % solve for fourier coeff (An)
c_L = x(1)*pi*AR; % coefficient of lift

delta = 0; % induced drag factor delta preallocation
for j=2:N % summing loop for odd indices
    delta = delta+(2*j-1)*(x(j)/x(1))^2; % sum
end
e = 1/(1+delta); % find e
c_Di = c_L^2/(pi*AR*e); % induced drag calculation
end

