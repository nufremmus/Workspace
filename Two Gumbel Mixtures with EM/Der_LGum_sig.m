%computes the derivative of log of Gumbel distribution wrt sigma
function y = Der_LGum_sig(mu,sig,x)

	y = -(x-mu)/(sig.^2) + exp((x-mu)/sig).*((x-mu)/(sig.^2)) - 1/sig;
		