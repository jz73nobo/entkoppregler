function [x] = norminv(p, mu, sigma, ~, ~)
	%NORMINV Inverse kumulierte Wahrscheinlichkeitsfunktion der Standardnormalverteilung N(mu, sigma^2)
	%	Input:
	%		p:		Wahrscheinlichkeit
	%		mu:		Erwatungswert
	%		sigma:	Varianz
	%		pcov:	?
	%		alpha:	?
	%	Output:
	%		x:		Wert für den p = P(x, N(mu, sigma^2)) gilt
	%#codegen
	if nargin > 3
		error('norminv', '%d Argumente werden nur mit der Statistic Toolbox unterstützt.', nargin);
	end
	x = emlnorminv(p, mu, sigma);
end