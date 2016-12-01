%%% In our simulation, we assume that one mixture is mainly one Gumbel and the other one mixture is mainly the other Gumbel distribution. The two Gumbel distributions are reasonably different.
% %testing parameters
mu1 = 2;
sig1 = 3;
n1 = 1000;

mu2 = 10;
sig2 = 4;
n2 = 1000;

alpha = 0.1;
beta = 0.9;
% %generate simulation data
U1 = rand(1,n1);
I1 = U1<alpha;
U2 = rand(1,n2);
I2 = U2<beta;

x1 = evrnd(mu1,sig1,1,n1);
x2 = evrnd(mu2,sig2,1,n1);
X1 = I1.*x1+(1-I1).*x2;

x1 = evrnd(mu1,sig1,1,n2);
x2 = evrnd(mu2,sig2,1,n2);
X2 = I2.*x1+(1-I2).*x2;

% true likelihood for two mixtures
True_ll = ((I1)*log(evpdf(X1,mu1,sig1))'+ (1-I1)*log(evpdf(X1,mu2,sig2))')/n1 + ((I2)*log(evpdf(X2,mu1,sig1))'+ (1-I2)*log(evpdf(X2,mu2,sig2))')/n2;

% fprintf ('True_ll is %f\n',True_ll);

[alpha_p,beta_p,mu1_p,sig1_p,mu2_p,sig2_p] = EM_Gumbel(X1,X2, 'moment', ' ');

%output
fprintf( 'alpha_p %f, beta_p %f, mu1_p %f, sig1_p %f, mu2_p %f, sig2_ %f \n' , alpha_p,beta_p,mu1_p,sig1_p,mu2_p,sig2_p);

%plot
binWidth = 0.05;
x = -40:binWidth:40;

%plot the simulated data
 figure
 hold on
plot(x,(alpha*evpdf(x,mu1,sig1)+(1-alpha)*evpdf(x,mu2,sig2)),'r','DisplayName','Simulated Mixture 1')
plot(x,(beta*evpdf(x,mu1,sig1)+(1-beta)*evpdf(x,mu2,sig2)),'b','DisplayName','Simulated Mixture 2');

%plot the EM result
plot(x,(alpha_p*evpdf(x,mu1_p,sig1_p)+(1-alpha_p)*evpdf(x,mu2_p,sig2_p)),'m--','DisplayName','Approx Mixture A')
plot(x,(beta_p*evpdf(x,mu1_p,sig1_p)+(1-beta_p)*evpdf(x,mu2_p,sig2_p)),'g--','DisplayName','Approx Mixture B');

legend('show')

saveas(gcf,'testing','pdf');