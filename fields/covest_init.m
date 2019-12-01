function [covf, par0, par_start] = covest_init(covf, par0, par_start, d_max, s2)
% covest_init Initial values for covest_ml and covest_ls
%
% Internal usage to set suitable initial values for covest_ml and covest_ls

%Determine covariance function
switch covf
	case 'cauchy'
		n_pars = 3;
		covf = @(d,x,I,J) cauchy_covariance(d,x(1),x(2),x(3),I,J);
		range_par = @(range, x) range/sqrt(10^(1/x)-1);
	case 'exponential'
		n_pars = 2;
		covf = @(d,x,I,J) exponential_covariance(d,x(1),x(2),I,J);
		range_par = @(range) 2/range;
	case 'gaussian'
		n_pars = 2;
		covf = @(d,x,I,J) gaussian_covariance(d,x(1),x(2),I,J);
		range_par = @(range) range;
	case 'matern'
		n_pars = 3;
		covf = @(d,x,I,J) matern_covariance(d,x(1),x(2),x(3),I,J);
		range_par = @(range, nu) sqrt(8*nu)/range;
	case 'spherical'
		n_pars = 2;
		covf = @(d,x,I,J) spherical_covariance(d,x(1),x(2),I,J);
		range_par = @(range) range/0.73;
	otherwise
		error('Unknown covariance function.')
end
%add nugget to number of parameters
n_pars = n_pars+1;
%check size of parameter vectors (or set suitable defaults)

%default for par0
if isempty(par0), par0=zeros(n_pars,1); end
%ensure that par0 is column vector
par0 = par0(:);
%check size of par0
if length(par0)~=n_pars, error('par0 should be of length %d', n_pars); end

%default for par_start
if isempty(par_start)
	par_start = par0;
	%rough estimate fo field variance
	par_start(1) = s2;
	%nugget should be one fifth of total variance
	par_start(end) = par_start(1)/5;
	%range should be about half of max distance
	if n_pars==3 %(exp, gaus, spher)
		par_start(2) = range_par(d_max/2);
	else %(matern, cauchy) first set defautl for nu,x
		if par_start(3)==0, par_start(3)=1; end
		par_start(2) = range_par(d_max/2, par_start(3));
	end
end
%ensure that par_start is column vector
par_start = par_start(:);
%check size of par_start
if length(par_start)~=n_pars
	error('par_start should be of length %d', n_pars);
end
