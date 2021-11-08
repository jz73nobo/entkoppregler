function [c_K, ceq_K, c_D, ceq_D, c_F, ceq_F, gradc_K, gradceq_K, gradc_D, gradceq_D, gradc_F, gradceq_F] = nonlin_constraint(K, D, F)
	c_K			= [];
	ceq_K		= [];
	c_D			= [];
	ceq_D		= [];
	c_F			= [];
	ceq_F		= [];
	gradc_K		= [];
	gradceq_K	= [];
	gradc_D		= [];
	gradceq_D	= [];
	gradc_F		= [];
	gradceq_F	= [];

	[p, n] = size(K);

	k = reshape(K, numel(K), 1);
	d = reshape(D, numel(D), 1);
	f = reshape(F, numel(F), 1);

	ceq_K = k(1) - k(2);
	gradceq_K = [1 1 zeros(1, numel(k) - 2)];
	gradceq_K = reshape(gradceq_K, p, n);
	% grad has to be p X n
end