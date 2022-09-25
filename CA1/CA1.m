clear all;clc;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CA01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CJ Kennedy - 109408903
% ASEN 3111 - Aerodynamics
%
% Contains Problem 1 and Problem 2
%
% Problem 1 needs functions simp.m and trap.m
% Problem 2 needs functions trap_question2.m
%
% Problem 1 produces two required plots and two plots used for error
% analysis. Also prints 3 lines
% Problem 2 prints 5 lines

%%
%%%%%% Problem 1 %%%%%%
%                     %
% assume: ideal flow (incompressible and inviscid) over rotating cyclinder
% given: Vinf, T, R, t (theta)
% find: 1-5
%                     %

% 1 - sectional coefficients of lift and drag
syms t R Vinf
T = -pi*R*Vinf;
Cp_fun = 1 - (4*sin(t)^2+2*T*sin(t)/(pi*R*Vinf)+(T/(2*pi*R*Vinf))^2);
cl_actual = -(1/2)*int(Cp_fun*sin(t),0,2*pi);
cd_actual = -(1/2)*int(Cp_fun*cos(t),0,2*pi);
fprintf(["Problem 1:"]); fprintf('\n')
fprintf(['Analytical sectional lift and drag coefficient are %.2f' ...
    ' and %.2f.'],cl_actual,cd_actual); fprintf('\n')
% 2 - Trapezoidal Rule
num = 10;
Cp_fun_cl = (1 - (4*sin(t)^2-2*sin(t)+(1/2)^2))*sin(t);
Cp_fun_cd = (1 - (4*sin(t)^2-2*sin(t)+(1/2)^2))*cos(t);
[cl_trap,cd_trap] = trap(Cp_fun_cl,Cp_fun_cd,num);
N = 1:num;
figure(1)
grid on
hold on
plot(N,cl_trap)
plot(N,cd_trap)
xlabel('Number of Panels (N)')
xticks(1:10)
ylabel('Sectional Coefficients')
legend('c_l(\theta)','c_d(\theta)')
title('Sectional Coefficients from Composite Trapezoidal Rule')
% 3 - Simpson's Rule
num = 10;
Cp_fun_cl = (1 - (4*sin(t)^2-2*sin(t)+(1/2)^2))*sin(t);
Cp_fun_cd = (1 - (4*sin(t)^2-2*sin(t)+(1/2)^2))*cos(t);
[cl_simp,cd_simp] = simp(Cp_fun_cl,Cp_fun_cd,num);
N = 1:num;
figure(2)
grid on
hold on
plot(N,cl_simp)
plot(N,cd_simp)
xlabel('Number of Panels (N)')
xticks(1:10)
ylabel('Sectional Coefficients')
legend('c_l(\theta)','c_d(\theta)')
title('Sectional Coefficients from Composite Simpsons Rule')
% 4 - Trap Error
trap_error = zeros(2,num);
for i=1:num
    trap_error(1,i) = -(cl_actual-cl_trap(i))*100;
    trap_error(2,i) = (cd_actual-cd_trap(i))*100;
end
figure(3)
hold on
plot(1:num,trap_error(1,:))
plot(1:num,trap_error(2,:))
xlabel('N')
ylabel('Error (in Percent)')
title('Trap Error vs N')
fprintf(['Trapezoidal: 1 percent relative error for sectional lift ' ...
    'coefficient requires N = 3 panels']); fprintf('\n')
% 5 - Simp Error
simp_error = zeros(2,num);
for i=1:num
    simp_error(1,i) = -(cl_actual-cl_simp(i))*100;
    simp_error(2,i) = -(cd_actual-cd_simp(i))*100;
end
figure(4)
hold on
plot(1:num,simp_error(1,:))
plot(1:num,simp_error(2,:))
xlabel('N')
ylabel('Error (in Percent)')
title('Simpsons Error vs N')
fprintf(['Simpsons: 1 percent relative error for sectional lift' ...
    ' coefficient requires N = 3 panels']); fprintf('\n')
fprintf('\n')
%%
%%%%%% Problem 2 %%%%%%
%                     %
% given: NACA 0012 airfoil @ 9deg angle of attack w
% find: L,D per unit span and error 
%                     %

% Given Data
c = 2;
alpha = 9;
Vinf = 60;
rho_inf = 1;
p_inf = 85.5*10^3;
q_inf = (1/2)*rho_inf*Vinf^2;
% Load Data
load Cp.mat
bigN = 5000; % num of iteration
[N,A] = trap_question2(Cp_upper,Cp_lower,bigN,q_inf,p_inf,c); % function
L = N*cosd(alpha)-A*sind(alpha);
D = N*sind(alpha)+A*cosd(alpha);
fprintf(["Problem 2:"]); fprintf('\n')
fprintf('Lift per unit span is %f',L); fprintf('\n')
fprintf('Drag per unit span is %f',D); fprintf('\n')


L = 3879.542070; D = 4.776080; % calculated with bigN = 100000
% Error Calculations n
n = 48;
[N,A] = trap_question2(Cp_upper,Cp_lower,n/2,q_inf,p_inf,c); 
L_prime = N*cosd(alpha)-A*sind(alpha);
error = (L-L_prime)/L*100;
fprintf(['10%% Relative Error ~ requires n = %.f points which' ...
    ' generates an error of %.3f%%'],n,error);
fprintf('\n')

n = 336;
[N,A] = trap_question2(Cp_upper,Cp_lower,n/2,q_inf,p_inf,c); 
L_prime = N*cosd(alpha)-A*sind(alpha);
error = (L-L_prime)/L*100;
fprintf(['1%% Relative Error ~ requires n = %.f points which' ...
    ' generates an error of %.3f%%'],n,error);
fprintf('\n')

n = 1950;
[N,A] = trap_question2(Cp_upper,Cp_lower,n/2,q_inf,p_inf,c); 
L_prime = N*cosd(alpha)-A*sind(alpha);
error = (L-L_prime)/L*100;
fprintf(['1/10%% Relative Error ~ requires n = %.f points which ' ...
    'generates an error of %.3f%%'],n,error);
fprintf('\n')
