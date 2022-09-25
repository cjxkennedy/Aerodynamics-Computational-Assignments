%% ASEN 3111 - Computational Assignment 3 - Main Script
% 
% Includes Solutions/Code for Three Problems:
% (1) - Comparison of the Lift Generated by Cambered Airfoils
% (2) - Analysis of Approximate Cessna 150 Wing Performance
% (3) - Validation of Span Efficiency Factor
%
% Where
% (1) Produces maximized subplot for NACA 0012, 2412, 4412 lift slopes
% (1) Prints lift slope, zero lift angle of attack for above airfiols
% (1) Included in above is corresponding thin airfoild theory solutions
% (2) Prints Lift and Induced Drag with given parameters
% (3) Produces Span Efficiency Factors vs Taper Ratio Plot four 4 AR's
%
% Author: CJ Kennedy
% Requires functions PLLT.m, NACA_Airoils, Vortex_Panel.m,
% aLoForThinCamber.m
% Date: 4/3/22
%
%% Housekeeping
clear all; clc; close all;
% tic
%% Problem #1
N = 70; % solution from CA2
c = 2; % wingspan
alpha = linspace(-15,15,N); % vary angle of attack (deg)
% NACA 0012
[x1,y1,aLo,dCLda] = NACA_Airfoils(0,0,12/100,c,N);
c_l = (alpha-aLo)*dCLda; % calculate c_l
aLo_00 = 0; % hardcoded zero lift angle of attack for sym thin airfoil
c_l_00 = (alpha-aLo_00)*dCLda;
set(gcf,'WindowState','maximized')
subplot(3,1,1)
hold on
plot(alpha,c_l)
plot(alpha,c_l_00) % plot thin airfoil theory 0012
xlabel('Angle of Attack (deg)'); ylabel('c_l');
title('c_l vs \alpha for NACA 0012')
legend('NACA 0012','Corresponding Thin-Airfoil Theory')
fprintf('Problem 1 -- Comparison of the Lift Generated by Cambered Airfoils. \r')
fprintf('For the NACA 0012, the lift slope = %.5f 1/deg and zero-lift %c = %.4f deg. \r',dCLda,(9082),aLo); 
fprintf('According to thin airfoil theory, lift slope = %.5f 1/deg and zero-lift %c = %.4f deg. \r',pi^2/90,(9082),aLo_00); 
% NACA 2412
[x2,y2,aLo,dCLda] = NACA_Airfoils(2/100,4/10,12/100,c,N);
c_l = (alpha-aLo)*dCLda;
aLo_24 = aLoForThinCamber(2/100,4/10,c);
c_l_24 = (alpha-aLo_24)*dCLda;
subplot(3,1,2)
hold on
plot(alpha,c_l)
plot(alpha,c_l_24) % plot thin airfoil theory 2412
xlabel('Angle of Attack (deg)'); ylabel('c_l');
title('c_l vs \alpha for NACA 2412')
legend('NACA 2412','Corresponding Thin-Airfoil Theory')
fprintf(['For the NACA 2412, the lift slope = %.5f 1/deg and ' ...
    'zero-lift %c = %.4f deg. \r'],dCLda,(9082),(aLo));
fprintf(['According to thin airfoil theory, lift slope = %.5f ' ...
    '1/deg and zero-lift %c = %.4f deg. \r'],pi^2/90,(9082),aLo_24); 
% NACA 4412
[x3,y3,aLo,dCLda] = NACA_Airfoils(4/100,4/10,12/100,c,N);
c_l = (alpha-aLo)*dCLda;
aLo_44 = aLoForThinCamber(4/100,4/10,c);
c_l_44 = (alpha-aLo_44)*dCLda;
subplot(3,1,3)
hold on
plot(alpha,c_l_44) % plot thin airfoil theory 4412
plot(alpha,c_l);
xlabel('Angle of Attack (deg)'); ylabel('c_l');
title('c_l vs \alpha for NACA 4412')
legend('NACA 4412','Corresponding Thin-Airfoil Theory')
fprintf(['For the NACA 4412, the lift slope = %.5f 1/deg and' ...
    ' zero-lift %c = %.4f deg. \r'],dCLda,(9082),(aLo));
fprintf(['According to thin airfoil theory, lift slope = %.5f ' ...
    '1/deg and zero-lift %c = %.4f deg. \r'],pi^2/90,(9082),aLo_44); 
sgtitle('Lift Slopes for Airfoils')
%% Problem #2
% Given Parameters
b = 33+4/12;  % span (ft)
c_t = 3+8.5/12; % chord @ tip (ft)
c_r = 5+4/12; % chord @ foot (ft)
geo_r = 1; % geo angle of attack @ root (deg)
geo_t = 0;  % geo angle of attack @ tip (deg)
Vinf = 42.184; % cruise speed (m/s)
[~,~,~, rho] = atmosisa(3048); % density at 10000ft / 3048 m
qinf = (1/2)*rho*Vinf^2; % dynamic pressure (SI Units)
% Convert to Metric (ft -> meters)
b = 0.3048*b;
c_t = 0.3048*c_t;
c_r = 0.3048*c_r;
S = b/2*(c_r+c_t); % planform area

N = 100; % spacing for NACA_Airfoils and VortexPanel
% Find aLo's and lift slope for airfoil at root and tip
[x1,y1,aero_r,a0_r] = NACA_Airfoils(2/100,4/10,12/100,c_r,N);
[x2,y2,aero_t,a0_t] = NACA_Airfoils(0,0,12/100,c_t,N); 
aero_t = 0; % set aLo @ tip = 0 to avoid small error
aero_r = deg2rad(aero_r); geo_r = deg2rad(geo_r); % convert to rad
a0_r = a0_r*180/pi; a0_t = a0_t*180/pi; % convert 1/deg to 1/rad
% Find e and coefficients of lift and induced drag w/ given parameters
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
% Find Lift and Drag and print
L = c_L*qinf*S;
D = c_Di*qinf*S;
fprintf(['\rProblem 2 -- Analysis of Approximate Cessna ' ...
    '150 Wing Performance. \r'])
fprintf('At given conditions, L = %.2f N and Di = %.2f N. \r', L, D);
fprintf(['Alternatively, L = %.2f lbs and Di = %.2f lbs. ' ...
    '\r'], L*0.22481, D*0.22481);

% Error 10%
L_10 = L*2; D_10 = D*2; % set initial value of L, D > actual L, D
ind = 1; % preallocate N value
while L_10 > L*(1.1) && D_10 > D*(1.1) % loop while error isn't accepted
    ind = ind+1; % increment N
    [~,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,ind);
    % calc new L and D
    L_10 = c_L*qinf*S;
    D_10 = c_Di*qinf*S;
end
fprintf(['Number of odd terms for solutions to be within 10%% ' ...
    'error is %.f. \r'],ind);
% Error 1%
L_1 = L*2; D_1 = D*2; % set initial value of L, D > actual L, D
ind = 1; % preallocate N value
while L_1 > L*(1.01) && D_1 > D*(1.01) 
    ind = ind+1; % increment
    [~,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,ind);
    % calc new L and D
    L_1 = c_L*qinf*S;
    D_1 = c_Di*qinf*S;
end
fprintf(['Number of odd terms for solutions to be within 1%% ' ...
    'error is %.f. \r'],ind);
% Error 1/10%
L_01 = L*2; D_01 = D*2; % set initial value of L, D > actual L, D
ind = 1; % preallocate N value
while L_01 > L*(1.001) && D_01 > D*(1.001)
    ind = ind+1;
    [~,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,ind);
    % calc new L and D
    L_01 = c_L*qinf*S;
    D_01 = c_Di*qinf*S;
end
fprintf(['Number of odd terms for solutions to be within 0.1%% ' ...
    'error is %.f. \r'],ind);
%% Problem #3
n = 50; % variable for spacing
AR=[4 6 8 10]; % varied aspect ratio
taper = linspace(0,1,n); % vary taper ratio 0 to 1
c_r = 1; % set constant root chord
c_t = taper*c_r; % vary tip chord w/ taper ratio constraints
figure()
for j = 1:4 % loop thru AR's
    b = AR(j)/2*(c_t+c_r); % calculate span
    e = zeros(1,n); % preallocate / reset e
    for i=1:n % loop through c_t / b to calc unique e's
        [e(i),~,~] = PLLT(b(i),2*pi,2*pi,c_t(i),c_r,0,0,1,1,n);
    end
    delta = 1./e-1; % induced drag factor calculation
    hold on
    plot(taper,e) % plot 
end
xlabel('c_t/c_r'); ylabel ('e')
legend('AR = 4' ,'AR = 6', 'AR = 8', 'AR = 10')
title('Span Efficiency Factor vs Taper Ratio')
%toc