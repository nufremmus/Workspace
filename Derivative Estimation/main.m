p = genpath(pwd);
addpath(p);
%------------------------Setup----------------------------------
k=3;

a = 0.25;
b = 1 ;
N = 500;
var = 0.01;

dx = (b-a)/N;
x=a:dx:b; x=x';
%------------------------start----------------------------------
h = @testfunction; %h is the function handle for true function of this model


y = distribution(a,b,N,h,var);%y is the data with noise

[ty,tdy] = testfunction(x); %ty and tdy are the true data and true derivative

G= Coeff_Matrix(k);
%Estimate derivative
dy = derivative_interior(a,b,N,y,k,'Penrose');

%plot
plot(x,dy,'or',x,tdy,'-g');
title(sprintf('plot_k=%d',k));
saveas(gca,sprintf('plot_k=%d',k),'epsc');
