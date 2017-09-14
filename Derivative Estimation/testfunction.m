function [ y,dy ]= testfunction(x)

y = sqrt((1-x).*x).*sin(2.1*pi./(x+0.05));
dy = (1-2*x)./(2*sqrt(x.*(1-x))).*sin(2.1*pi./(x+0.05)) - 2.1*pi./((x+0.05).^2).*sqrt((1-x).*x).*cos(2.1*pi./(x+0.05)); %dy = -inf at x =1

% y = sin(x);
% dy = cos(x);