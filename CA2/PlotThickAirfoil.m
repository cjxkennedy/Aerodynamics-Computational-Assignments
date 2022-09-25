function [] = PlotThickAirfoil(c,alpha,V_inf,p_inf,rho_inf,N)

%PlotThickAirfoil func plots (1) stream lines, (2) equipotential lines, 
% and (3) pressure contours
%   Approximations for flow about a thick and symmetric airfoil. Inputs are
%   physical parameters and N - discrete vortices used in vortex sheet
%   approximation. 
%
%   Heavily adapted from PlotThinAirfoil.m and Lifting_Cylinder.m
%
%   Requires Vortex_Panel.m 
%
%   Author: CJ Kennedy
%   Date: 2/27/21

%% Define Grid Points, Domain, and Mesh
nx=500; % steps in the x direction
ny=500; % steps in the y direction
xmin=-1;
xmax=2.5;
ymin=-1;
ymax=1;
[x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));
%% Setup Vectors and Call Vortex Panel
A = linspace(c,0,N/2); % first half of x_b
B = flip(A(1:end-1)); % second half of x_b
C = 12/100*(c/.2)*(.2969*sqrt(A./c)-.126*(A./c)-.3516*(A./c).^2+ ...
   .2843*(A./c).^3-.1036*(A./c).^4); % first half of y_b 
D = 12/100*(c/.2)*(.2969*sqrt(B./c)-.126*(B./c)-.3516*(B./c).^2+ ...
   .2843*(B./c).^3-.1036*(B./c).^4); % second half of y_b 
x_b = [A B]; % catinate x vectors
y_b = [-C D]; % catinate y vectors (note -C to flip sign)
[c_l, c_p, circ, x_i, y_i] = Vortex_Panel(x_b,y_b,V_inf,rad2deg(alpha),0);
%% Setup For-Loop Variables
% Preallocate Arrays
Phi = zeros(nx,ny);
Psi = zeros(nx,ny);
V1 = zeros(nx,ny);
V2 = zeros(nx,ny);
%% Loop through vortices N 
for i=1:N-2 % 
    T = circ(i); % circulation for current vortex
    x1 = x_i(i); % x position
    y1 = y_i(i); % y position
    t = mod(atan2(y-y1,x-x1),2*pi); % theta (polar coords)
    r = sqrt((x-x1).^2+(y-y1).^2); % radius (polar coords)
    % Streamline and Velocity Potential Calculation 
    Phi = Phi-T/(2*pi)*t; 
    Psi = Psi+T/(2*pi)*log(r);
    % Velocity Calculations
    Vx = (T*sin(t))./(2*r*pi);
    Vy = -(T*cos(t))./(2*r*pi);
    V1 = V1+Vx;
    V2 = V2+Vy;
end
% Calculate psi and phi for uniform stream
Psi_U0 = V_inf*y*cos(alpha)-V_inf*x*sin(alpha);
Phi_U0 = V_inf*x*cos(alpha)+V_inf*y*sin(alpha);
% Also calculate x and y velocities from uniform stream
Ux = V_inf*cos(alpha);
Uy = V_inf*sin(alpha);
% Find total streamfunction and velocity potential
PHI_T = Phi_U0+Phi;
PSI_T = Psi_U0+Psi;
%% Plotting equipotential lines
figure() % setup subplot that's maximized and has title
set(gcf,'WindowState','maximized')
sgtitle('Problem 2')
subplot(2,2,1)
hold on
title('Equipotential Lines')
levels = linspace(-100,150,60)';
contour(x,y,PHI_T,levels,'LineWidth',1.5)
axis equal
axis([xmin xmax ymin ymax]) % set axis bounds for plots
fill(x_i,y_i,'w') % shade airfoil region
plot(x_i,y_i,'Color','Black','LineWidth',1.5) % outline airfoil region
ylabel('y [m]'); xlabel('x [m]')
%% Plotting stream lines
% Determine color levels for contours
levmin = PSI_T(1,nx); % defines the color levels
levmax = PSI_T(ny,nx/4);
levels = linspace(levmin,levmax,50)';
% Plot streamfunction at levels
subplot(2,2,2)
hold on
title('Streamlines')
contour(x,y,PSI_T,levels,'LineWidth',1.5)
axis equal
axis([xmin xmax ymin ymax])  % set axis bounds for plots
plot(x_i,y_i,'Color','Black','LineWidth',1.5)
ylabel('y [m]'); xlabel('x [m]')
%% Plot pressure contours
V1 = V1 + Ux; % add uniform vel x comp with x components of vortex's
V2 = V2 + Uy; % add uniform vel y comp with y components of vortex's
V = sqrt(V1.^2+V2.^2); % find magnitude
Cp = 1 - (V/V_inf).^2; % compute coefficient of pressure (Ber's eq)
subplot(2,2,[3,4])
hold on
title('Pressure Countours')
levels = linspace(-5,1,30)';
contourf(x,y,Cp,levels,'LineColor','None')
axis equal
axis([xmin xmax ymin ymax]) % set axis bounds for plots
fill(x_i,y_i,'white') % shade airfoil region
plot(x_i,y_i,'Color','Black','LineWidth',1.5) % outline airfoil region
ylabel('y [m]'); xlabel('x [m]')
colorbar