function [wNew,PNew] = exrnq_update( ...
    w,P,uvec,d,A,beta,q,j)

% ============================================================
% EX-RNQ
% One iteration of the Extended Recursive
% Nonquadratic algorithm
% ============================================================

L = length(w);

%% Nonquadratic coefficient

c = (2*j - 1)*d^(2*j - 2);

%% Numerical protection when d is close to zero

cMin = 1e-12;

if abs(c) < cMin
    c = cMin;
end

invc = 1/c;

%% R quantity

R = invc ...
    + (1/beta)*(uvec.'*P*uvec);

%% Gain vector

G = ((1/beta)*A*P*uvec)/R;

%% Modified innovation

epsilon = ...
    d/(2*j - 1) ...
    - uvec.'*w;

%% Weight update

wNew = A*w + G*epsilon;

%% P update

PNew = ...
    (1/beta)*A*P*A.' ...
    - (1/beta)*G*G.'*R ...
    + invc*q*eye(L);

%% Numerical symmetry

PNew = (PNew + PNew.')/2;

end