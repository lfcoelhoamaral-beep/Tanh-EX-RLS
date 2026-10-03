function [wNew,VNew,state,K,qRLM,e] = rlm_update( ...
    w,V,uvec,d,state,lambda,lambdaSigma,Nw)

% ============================================================
% RLM_UPDATE
%
% Recursive Least M-Estimate (RLM) algorithm
% Zou, Chan and Ng, IEEE Signal Processing Letters, 2000
%
% Hampel three-part redescending M-estimate
% ============================================================

    L = length(w);

    %% ========================================================
    % A priori estimation error
    % =========================================================

    e = d - w.'*uvec;

    %% ========================================================
    % Update window of squared errors
    % A_e(n) = {e^2(n), ..., e^2(n-Nw+1)}
    % =========================================================

    e2 = e^2;

    if isempty(state.err2Buffer)

        state.err2Buffer = e2;

    else

        state.err2Buffer = ...
            [e2; state.err2Buffer(:)];

        if length(state.err2Buffer) > Nw
            state.err2Buffer = ...
                state.err2Buffer(1:Nw);
        end

    end

    %% ========================================================
    % Robust recursive variance estimate
    % =========================================================

    C1 = 1.483*(1 + 5/(Nw-1));

    medE2 = median(state.err2Buffer);

    if ~state.initialized

        % The original paper initializes sigma_hat^2(0)
        % from the squared desired signal.
        state.sigma2 = max(d^2,1e-12);
        state.initialized = true;

    end

    sigma2New = ...
        lambdaSigma*state.sigma2 ...
        + C1*(1-lambdaSigma)*medE2;

    sigma2New = max(sigma2New,1e-12);

    sigmaHat = sqrt(sigma2New);

    %% ========================================================
    % Hampel thresholds
    % =========================================================

    xi     = 1.960*sigmaHat;
    Delta1 = 2.240*sigmaHat;
    Delta2 = 2.576*sigmaHat;

    absE = abs(e);

    %% ========================================================
    % Robust weighting function
    %
    % q(e) = psi(e)/e
    % =========================================================

    if absE < xi

        qRLM = 1;

    elseif absE < Delta1

        qRLM = xi/absE;

    elseif absE < Delta2

        qRLM = ...
            xi*(Delta2-absE) / ...
            (absE*(Delta2-Delta1));

    else

        qRLM = 0;

    end

    %% ========================================================
    % M-estimate gain vector
    %
    % K(n) =
    % q(e)V(n-1)X(n) /
    % [lambda + q(e)X^T(n)V(n-1)X(n)]
    % =========================================================

    denominator = ...
        lambda ...
        + qRLM*(uvec.'*V*uvec);

    K = ...
        qRLM*V*uvec / denominator;

    %% ========================================================
    % Weight update
    % =========================================================

    wNew = ...
        w + K*e;

    %% ========================================================
    % Inverse correlation matrix update
    %
    % V(n) = lambda^{-1}
    %        [I - K(n)X^T(n)]V(n-1)
    % =========================================================

    VNew = ...
        (1/lambda) * ...
        (eye(L)-K*uvec.')*V;

    %% Numerical symmetry

    VNew = (VNew + VNew.')/2;

    %% Update state

    state.sigma2 = sigma2New;

end