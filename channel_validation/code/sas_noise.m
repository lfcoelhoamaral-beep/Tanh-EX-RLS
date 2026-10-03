function x = sas_noise(alpha,sigma,N)
% SAS_NOISE Generates symmetric alpha-stable noise.
%
% Characteristic function:
% phi(t) = exp(-sigma^alpha * |t|^alpha)
%
% Inputs:
%   alpha : characteristic exponent, 0 < alpha <= 2
%   sigma : scale parameter, sigma > 0
%   N     : number of samples
%
% Output:
%   x     : N-by-1 symmetric alpha-stable random vector

if alpha <= 0 || alpha > 2
    error('alpha must satisfy 0 < alpha <= 2.');
end

if sigma <= 0
    error('sigma must be positive.');
end

if N <= 0 || floor(N) ~= N
    error('N must be a positive integer.');
end

% Gaussian limiting case: alpha = 2
if abs(alpha - 2) < 1e-12
    x = sqrt(2)*sigma*randn(N,1);
    return;
end

% Auxiliary random variables
V = pi*(rand(N,1) - 0.5);
W = -log(rand(N,1));

% Cauchy case: alpha = 1
if abs(alpha - 1) < 1e-12
    x = sigma*tan(V);
    return;
end

% Chambers-Mallows-Stuck method
term1 = sin(alpha*V) ./ ...
    (cos(V)).^(1/alpha);

term2 = ...
    (cos((1-alpha)*V) ./ W) ...
    .^((1-alpha)/alpha);

x = sigma .* term1 .* term2;

end