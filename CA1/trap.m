function [cl,cd] = trap(func1,func2,num)
% function for using composite trapezoidal rule
% inputs are the functions required and # of panels to iterate over (num)
% outputs are the two appropiate approximations
syms t
for i=1:num
    N = i;
    start = 0; stop = 2*pi;
    data = linspace(start,stop,N+1);
    approx1 = 0;
    approx2 = 0;
    for k=1:N
        approx1 = approx1+(data(k+1)-data(k))*(subs(func1,t,data(k+1))+subs(func1,t,data(k)))/2;
        approx2 = approx2+(data(k+1)-data(k))*(subs(func2,t,data(k+1))+subs(func2,t,data(k)))/2;
    end
    cl(i) = -(1/2)*approx1;
    cd(i) = -(1/2)*approx2;
end
end

