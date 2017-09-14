p = genpath(pwd);
addpath(p);
%------------------------Setup----------------------------------

load Gaussian_mixProp05_loc2
k = 5;  %parameter for generating data in knnsearch
L = 12;  %L is number of symmetric intervals used;L could be any integer >=2
dx=0.01;  %the space between x; in your case is alpha
a=0; b=1; %the range of x is [a,b]

x = a:dx:b;
N =(b-a)/dx;
%-------------------generating equi-spaced data---------------
alphas1 = x;
IDX=knnsearch(alphas',alphas1','K',k);
fs_est=median(fs(IDX),2)';

y = fs_est;
%------------------------start----------------------------------
dy = zeros(1,N+1);

% compute derivative for the interior points; the boundary derivatives are not estimated
for j=L+1:N-L
   for i=1:L
	   coeff = 6*(i^2)/(L*(L+1)*(2*L+1));
	   fin_diff = (y(j+i)-y(j-i))/(2*i*dx); 
	   dy(j) = dy(j)+coeff*fin_diff;
	end
end
%plot
plot(x,dy,'or');
title(sprintf('Kris_Derivative_plot_k=%d_L=%d',k,L));
saveas(gca,sprintf('Kris_Derivative_plot_k=%d_L=%d',k,L),'epsc');
