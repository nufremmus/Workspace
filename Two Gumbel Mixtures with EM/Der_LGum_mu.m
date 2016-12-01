%computes the derivative of log of Gumbel distribution wrt mu
function y = Der_LGum_sig(mu,sig,x)

	y = -1/sig + 1/sig*exp((x-mu)./sig);
		