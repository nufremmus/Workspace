function y = distribution(a,b,N,h,var)

% var is the variance of the error ~ N(0,var)

dx = (b-a)/N;
x = a:dx:b; x=x';
[y,dy] = h(x); % the true model without noise

%generate the errors
err = zeros(N+1,1);
for i=1:N+1
	err(i) = normrnd(0,sqrt(var));
end

y = y+err;
