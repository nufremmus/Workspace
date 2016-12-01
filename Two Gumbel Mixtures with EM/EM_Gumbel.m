% We are given the following scenario
% g=alpha*f1+(1-alpha)*f2
% h=beta*f1+(1-beta)*f2
% We are given two sets of data and one is from g with n1 data and the other is from h with n2 data
% f1 and f2 are Gumbel distributions
% Now we want to use EM to find the proportions as well as the distributions

% The following is the outline of this algorithm
% 1. Initialize the classifiers for each x, e.g. K-mean
% repeat the following:
% 2. For 1st dataset apply the Expectation step
% 3. For 1st dataset apply M-step
% 4. For 2nd dataset apply E-step
% 5. For 2nd dataset apply M-step

% Alternative outline of this algorithm
% 1. Initialization of the membership
% repeat:
% 2. Apply E-step for two datasets independently
% 3. Apply M-step to the combined dataset


% References:
% SAR IMAGES AS MIXTURES OF GAUSSIAN MIXTURES,Peter Orbanz and Joachim M. Buhmann,Institute of Computational Science, ETH Zurich

% Note : our gumbel distributions follow the definition in matlab, that is the gumbel distribution have diminishing tails on the left hand side

function [alpha,beta,mu1,sig1,mu2,sig2] = EM_Gumbel(X1,X2, init_method, opt_method)

	%x1 amd x2 are row vectors
	
	n1 = length(X1);
	n2 = length(X2);
	
	%initialization of em parameters
	tolerance = 10^(-4);
	max_gd = 100; %max round for gradient descend
	learn_rate = 0.01;
	prev_val = 0; %previous em likelihood
	curr_val = 0; %current em likelihood	

	%initialization of em parameters
	tolerance_em = 10^(-4);
	max_em = 10000; % max round for EM algorithm
	prev_val_em = 0; %previous likelihood
	curr_val_em = 0; %current likelihood	

	%initialization
	switch init_method
		case 'random'
		 		rng(1);
			mu1 = rand;
			sig1 = 1+rand;
			mu2 = rand;
			sig2 = 1+rand;
			alpha = rand;
			beta = rand;
			
		case 'determined'
		
			mu1 = 0.5;
			sig1 = 1;
			mu2 = 3;
			sig2 = 2;
			%mixing value
			alpha = 0.1;
			beta = 0.9;
			
		case 'moment'
			rng;
			sig1 = sqrt(var(X1)*6/pi^2)+rand;
			mu1 = mean(X1) + 0.577216*sig1+rand;
			% mu1 = mode(X1);
			sig2 = sqrt(var(X2)*6/pi^2)+rand;
			mu2 = mean(X2) + 0.577216*sig2+rand;
			% mu2 = mode(X2);
			alpha = rand;
			beta = rand;

		otherwise
			disp('Not defined')
	end
	
	% EM algorithm
	for i = 1:max_em
		fprintf('round %d\n' , i);
		% E-step
		%compute membership probabilities for each data in X1 belonging in f1 and f2 respectively
		
		z1_X1 = alpha*evpdf(X1,mu1,sig1)./(alpha*evpdf(X1,mu1,sig1) + (1-alpha)*evpdf(X1,mu2,sig2)); %
		z2_X1 = 1-z1_X1;
		%Compute membership probabilities for each data in X2 belonging in f1 and f2 respectively
		z1_X2 = beta*evpdf(X2,mu1,sig1)./(beta*evpdf(X2,mu1,sig1) + (1-beta)*evpdf(X2,mu2,sig2)); %
		z2_X2 = 1-z1_X2;
		
		% M-step
		% Maximize the combined log likelihood of two mixtures which is
		% [(z1_X1).*log(evpdf(X1,mu1,sig1))+ (z2_X1).*log(evpdf(X1,mu2,sig2))]+[(z1_X2).*log(evpdf(X2,mu1,sig1))+ (z2_X2).*log(evpdf(X2,mu2,sig2))] 
		
		prev_val_em = (z1_X1*log(evpdf(X1,mu1,sig1))'+ z2_X1*log(evpdf(X1,mu2,sig2))')/n1 + ((z1_X2)*log(evpdf(X2,mu1,sig1))'+ (z2_X2)*log(evpdf(X2,mu2,sig2))')/n2;
		prev_val = prev_val_em;

		%% start gradient descend
		for j=1:max_gd
			 %fprintf('round (%d, %d)\n' , i,j);
		
			% the derivative of the above likelihood is
			% wrt sig1: 
			step_sig1 = z1_X1*(Der_LGum_sig(mu1,sig1,X1))'+(z1_X2)*(Der_LGum_sig(mu1,sig1,X2))';
			%wrt mu1:
			step_mu1 = z1_X1*(Der_LGum_mu(mu1,sig1,X1))'+(z1_X2)*(Der_LGum_mu(mu1,sig1,X2))';
			% wrt sig2: 
			step_sig2 = z2_X1*(Der_LGum_sig(mu2,sig2,X1))'+(z2_X2)*(Der_LGum_sig(mu2,sig2,X2))';
			%wrt mu2:
			step_mu2 = z2_X1*(Der_LGum_mu(mu2,sig2,X1))'+(z2_X2)*(Der_LGum_mu(mu2,sig2,X2))';
			
			%update parameters
			% mu1 = mu1 + learn_rate*step_mu1/norm(step_mu1);
			% sig1 = sig1 + learn_rate*step_sig1/norm(step_sig1);
			% mu2 = mu2 + learn_rate*step_mu2/norm(step_mu2);
			% sig2 = sig2 + learn_rate*step_sig2/norm(step_sig2);

			step = norm([step_mu1 step_mu2 step_sig1 step_sig2]);
			mu1 = mu1 + learn_rate*step_mu1/step;
			sig1 = sig1 + learn_rate*step_sig1/step;
			mu2 = mu2 + learn_rate*step_mu2/step;
			sig2 = sig2 + learn_rate*step_sig2/step;
			
			%get current likelihood
			curr_val = ((z1_X1)*log(evpdf(X1,mu1,sig1))'+ (z2_X1)*log(evpdf(X1,mu2,sig2))')/n1 + ((z1_X2)*log(evpdf(X2,mu1,sig1))'+ (z2_X2)*log(evpdf(X2,mu2,sig2))')/n2;

			%compare likelihood difference
			if curr_val - prev_val< tolerance
				break;
			end
			
			prev_val = curr_val;
			
		end
		
		%update the curr_val_em
		curr_val_em = curr_val;
		
		%output current LL
		% fprintf('Current LL is %f\n', curr_val);
		fprintf('Current parameters are: LL = %.5f mu1=%2.3f sig2 = %.3f mu2 = %2.3f sig2 = %.3f alpha=%.2f beta = %.2f \n', curr_val,mu1,sig1,mu2,sig2,alpha,beta);
		
		%check is current loglikelihood is infinity
		if(isinf(curr_val))
			fprintf('%f %f %f %f\n',(z1_X1)*log(evpdf(X1,mu1,sig1))',(z2_X1)*log(evpdf(X1,mu2,sig2))',(z1_X2)*log(evpdf(X2,mu1,sig1))',(z2_X2)*log(evpdf(X2,mu2,sig2))');
			fprintf('%f %f  %f\n',sum(z1_X2),sum(log(evpdf(X2,mu1,sig1))), sum(find((X2-mu1)/sig1 < -5)));
			break;
		end
		
		%re-initialize parameters
		if(isnan(curr_val))
			fprintf('\n!!!REINITIALIZING!!!\n', curr_val);
							
			switch init_method
				case 'random'
					rng(i)
					mu1 = rand;
					sig1 = 1+rand;
					mu2 = rand;
					sig2 = 1+rand;
					alpha = rand;
					beta = rand;
					
				case 'determined'
					mu1 = 0.5;
					sig1 = 1;
					mu2 = 3;
					sig2 = 2;
					%mixing value
					alpha = 0.5;
					beta = 0.5;
					
				case 'moment'
					rng(i)
					sig1 = sqrt(var(X1)*6/pi^2)+rand;
					mu1 = mean(X1) + 0.577216*sig1+rand;
					% mu1 = mode(X1);
					sig2 = sqrt(var(X2)*6/pi^2)+rand;
					mu2 = mean(X2) + 0.577216*sig2+rand;
					% mu2 = mode(X2);
					alpha = rand;
					beta = rand;
					
				case 'kmean'
					
				otherwise
					disp('Not defined')
			end
			
			continue;
		end
		
		%compare the EM-likelihood difference
		if curr_val_em - prev_val_em< tolerance_em
			fprintf('final loglikelihood is %f\n',curr_val_em);
			break;
		end
		
		% update the mixing values for each mixtures
		alpha = sum(z1_X1)/n1;
		beta = sum(z1_X2)/n2;
		
	end
end