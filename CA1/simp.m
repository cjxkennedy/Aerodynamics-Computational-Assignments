function [cl,cd] = trap(func1,func2,num)
% function for using composite Simpson's rule
% inputs are the functions required and # of panels to iterate over (num)
% outputs are the two appropiate approximations
syms t
for i=1:num
    N = i;
    start = 0; stop = 2*pi;
    data = linspace(start,stop,2*N+1);
    approx1 = 0;
    approx2 = 0;
    for k=1:N
        approx1 = approx1+(subs(func1,t,data(2*k+1))+4*subs(func1,t,data(2*k))+subs(func1,t,data(2*k+1)));
        approx2 = approx2+(subs(func2,t,data(2*k+1))+4*subs(func2,t,data(2*k))+subs(func2,t,data(2*k+1)));
    end
    h = (stop-start)/(2*N);
    cl(i) = -(1/2)*approx1*h/3;
    cd(i) = -(1/2)*approx2*h/3;
end
end

