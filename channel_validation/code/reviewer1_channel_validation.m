function reviewer1_channel_validation(stage)
% Comparable bounded-budget tuning for nonlinear channel output prediction.
% Usage: reviewer1_channel_validation('smoke'), ('validate'), then ('test').
% Requires sas_noise, exrnq_update and rlm_update on the MATLAB path.
% No global-optimality or physical-channel recovery claim.
if nargin == 0, stage = 'smoke'; end
stage = lower(char(stage));
assert(any(strcmp(stage,{'smoke','validate','test'})), 'Unknown stage.');
for helper = {'sas_noise','exrnq_update','rlm_update'}
    assert(exist(helper{1},'file')~=0, 'Missing helper: %s', helper{1});
end
root = fileparts(mfilename('fullpath'));
protocol.L=5; protocol.N=2000; protocol.change=1001;
protocol.windows={901:1000,1001:1200,1901:2000};
protocol.A=0.9999999368*eye(protocol.L);
protocol.Ts=0.8e-6; protocol.fD1=100; protocol.fD2=50;
protocol.methods={'EX-RLS','Tanh-EX-RLS','EX-RNQ','RLM'};
protocol.scenarios=[2 .1;1.5 .1;1.2 .1;1 .1];
protocol.K=24; protocol.validationMC=30; protocol.testMC=200;
protocol.validationSeedBase=410000; protocol.testSeedBase=1010000;
protocol.candidateSeed=83017;
protocol.scoreDefinition=['Equal mean across scenarios of the equal mean ' ...
    'of pre-change, transient and post-change MedSE(dB). Lower is better.'];
if strcmp(stage,'smoke')
    protocol.K=3; protocol.validationMC=3;
    protocol.scenarios=[2 .1;1.5 .1];
    out=fullfile(root,'results_channel_smoke');
else
    out=fullfile(root,'results_channel_full');
end
if ~exist(out,'dir'), mkdir(out); end
% Preserve caller RNG state. Generating candidates never changes data seeds.
previousRNG=rng; cleanup=onCleanup(@()rng(previousRNG)); %#ok<NASGU>
if strcmp(stage,'test')
    selectionFile=fullfile(out,'selection.mat');
    assert(exist(selectionFile,'file')==2, 'Run validate first.');
    loaded=load(selectionFile,'protocol','selected','bestIndex');
    protocol=loaded.protocol; selected=loaded.selected;
    bestIndex=loaded.bestIndex;
    assert(protocol.K==24 && protocol.validationMC==30, ...
        'Testing requires the full protocol, not smoke results.');
    fprintf('Independent test: parameters frozen; %d runs/scenario.\n',protocol.testMC);
    bank=make_bank(protocol,protocol.testMC,protocol.testSeedBase);
    rows=cell(0,11);
    for m=1:4
        fprintf('Testing %s...\n',protocol.methods{m});
        [metrics,curves,failed,issue]=evaluate(protocol,bank,m,selected{m});
        fprintf('Completed %s; numerical failures: %d\n',protocol.methods{m},sum(failed));
        for a=1:size(protocol.scenarios,1)
            rows(end+1,:)={protocol.methods{m},bestIndex(m), ...
                protocol.scenarios(a,1),protocol.scenarios(a,2),protocol.testMC, ...
                metrics(a,1),metrics(a,2),metrics(a,3),failed(a), ...
                recovery(curves(a,:),metrics(a,3),protocol.change),issue{a}}; %#ok<AGROW>
        end
        save(fullfile(out,sprintf('test_method_%d.mat',m)), ...
             'protocol','selected','bestIndex','metrics','curves','failed','issue');
    end
    summary=cell2table(rows,'VariableNames',{'Method','Candidate','Alpha','Sigma', ...
        'Runs','PreMedSE_dB','TransientMedSE_dB','PostMedSE_dB', ...
        'NumericalFailureCount','RecoverySamples','FirstFailure'});
    writetable(summary,fullfile(out,'independent_test_summary.csv'));
    save(fullfile(out,'independent_test_summary.mat'),'summary','protocol','selected');
    disp(summary); return
end
candidates=make_candidates(protocol);
save(fullfile(out,'protocol_and_candidates.mat'),'protocol','candidates');
write_candidates(protocol,candidates,fullfile(out,'candidate_parameters.csv'));
bank=make_bank(protocol,protocol.validationMC,protocol.validationSeedBase);
scores=inf(4,protocol.K); numericalFailures=zeros(4,protocol.K);
selected=cell(4,1); bestIndex=zeros(4,1); rows=cell(0,10);
for m=1:4
    for k=1:protocol.K
        [metrics,~,failed,issue]=evaluate(protocol,bank,m,candidates{m,k});
        numericalFailures(m,k)=sum(failed);
        % Any nonfinite trajectory disqualifies the candidate, not omitted.
        if all(failed==0) && all(isfinite(metrics(:)))
            scores(m,k)=mean(metrics(:));
        end
        for a=1:size(protocol.scenarios,1)
            rows(end+1,:)={protocol.methods{m},k,protocol.scenarios(a,1), ...
                protocol.scenarios(a,2),metrics(a,1),metrics(a,2),metrics(a,3), ...
                failed(a),issue{a},scores(m,k)}; %#ok<AGROW>
        end
        fprintf('%s candidate %d/%d: score %.4f, failures %d\n', ...
            protocol.methods{m},k,protocol.K,scores(m,k),numericalFailures(m,k));
        save(fullfile(out,'validation_checkpoint.mat'), ...
            'protocol','candidates','scores','numericalFailures','rows');
    end
    [bestScore,bestIndex(m)]=min(scores(m,:));
    assert(isfinite(bestScore),'All candidates failed for %s.',protocol.methods{m});
    % Exact ties select first index deterministically.
    selected{m}=candidates{m,bestIndex(m)};
end
validation=cell2table(rows,'VariableNames',{'Method','Candidate','Alpha','Sigma', ...
    'PreMedSE_dB','TransientMedSE_dB','PostMedSE_dB', ...
    'NumericalFailureCount','FirstFailure','AggregateScore_dB'});
writetable(validation,fullfile(out,'validation_scores.csv'));
save(fullfile(out,'selection.mat'),'protocol','selected','bestIndex','scores', ...
     'numericalFailures','validation');
for m=1:4
    fprintf('\nSelected %s: candidate %d, score %.4f\n', ...
        protocol.methods{m},bestIndex(m),scores(m,bestIndex(m)));
    disp(selected{m});
end
if strcmp(stage,'smoke')
    fprintf('\nSMOKE ONLY: do not use these selections as paper results.\n');
else
    fprintf('\nValidation complete. Run reviewer1_channel_validation(''test'').\n');
end
end

function C=make_candidates(p)
rng(p.candidateSeed,'twister'); C=cell(4,p.K);
base=struct('q',3.26e-7,'r0',1e-2,'initialP',1e3,'gamma',4, ...
    'beta',.995,'j',2,'lambda',.99,'lambdaSigma',.99,'Nw',14);
for m=1:4, C{m,1}=base; end
C{3,1}.q=3.26e-7; C{4,1}.initialP=20;
% Equal candidate budget. Same nuisance candidates for EX-RLS and Tanh.
% Random search in explicitly bounded ranges; baseline is always included.
for k=2:p.K
    c=base; c.q=10^(-8+4*rand); c.r0=10^(-3+2*rand);
    c.initialP=10^(1+3*rand);
    C{1,k}=c; C{2,k}=c;
    C{2,k}.gamma=10^(log10(.25)+rand*log10(16/.25));
    c=base; c.beta=.95+(.9995-.95)*rand;
    c.q=10^(-8+4*rand); c.initialP=10^(1+3*rand);
    orders=[2 3]; c.j=orders(randi(numel(orders))); C{3,k}=c;
    c=base; c.lambda=.95+(.9995-.95)*rand;
    c.lambdaSigma=.95+(.9995-.95)*rand;
    windows=[8 14 20 32]; c.Nw=windows(randi(numel(windows)));
    c.initialP=10^(1+3*rand); C{4,k}=c;
end
end

function write_candidates(p,C,path)
rows=cell(4*p.K,12); n=0;
for m=1:4
    for k=1:p.K
        c=C{m,k};n=n+1;
        % NaN denotes a parameter unused by that method.
        q=NaN;r0=NaN;g=NaN;b=NaN;j=NaN;la=NaN;ls=NaN;nw=NaN;
        if m<=3, q=c.q; end
        if m<=2, r0=c.r0; end
        if m==2, g=c.gamma; end
        if m==3, b=c.beta;j=c.j; end
        if m==4, la=c.lambda;ls=c.lambdaSigma;nw=c.Nw; end
        rows(n,:)={p.methods{m},k,k==1,c.initialP,q,r0,g,b,j,la,ls,nw};
    end
end
T=cell2table(rows,'VariableNames',{'Method','Candidate','SubmittedBaseline', ...
    'InitialP','q','r0','gamma','beta','j','lambda','lambdaSigma','Nw'});
writetable(T,path);
end

function bank=make_bank(p,MC,seedBase)
S=size(p.scenarios,1);bank=cell(S,MC);
for a=1:S
    for run=1:MC
        rng(seedBase+10000*a+run,'twister');
        u=randn(p.N+p.L-1,1);
        H=zeros(p.L,p.N);
        h=(randn(p.L,1)+1i*randn(p.L,1))/sqrt(2);
        for i=1:p.N
            if i<p.change,fD=p.fD1;else,fD=p.fD2;end
            rho=besselj(0,2*pi*fD*p.Ts);
            eta=(randn(p.L,1)+1i*randn(p.L,1))/sqrt(2);
            h=rho*h+sqrt(max(1-rho^2,0))*eta;
            H(:,i)=real(h);
        end
        v=sas_noise(p.scenarios(a,1),p.scenarios(a,2),p.N);
        U=zeros(p.L,p.N);target=zeros(1,p.N);
        for i=1:p.N
            idx=i+p.L-1;U(:,i)=u(idx:-1:idx-p.L+1);
            target(i)=tanh(U(:,i).'*H(:,i)+v(i));
        end
        bank{a,run}=struct('U',U,'d',target);
    end
end
end

function [metrics,curves,failed,issue]=evaluate(p,bank,method,c)
S=size(bank,1);MC=size(bank,2);metrics=nan(S,3);
curves=nan(S,p.N);failed=zeros(S,1);issue=repmat({''},S,1);
I=eye(p.L);A=p.A;
for a=1:S
    deviation=nan(MC,p.N);
    for run=1:MC
        w=zeros(p.L,1);P=c.initialP*I;
        state=struct('sigma2',0,'err2Buffer',[],'initialized',false);
        D=bank{a,run};
        for i=1:p.N
            u=D.U(:,i);d=D.d(i);
            prediction=w.'*u; % A priori output, BEFORE any coefficient update.
            try
                if method==4
                    [w,P,state]=rlm_update(w,P,u,d,state,c.lambda,c.lambdaSigma,c.Nw);
                elseif method==3
                    [w,P]=exrnq_update(w,P,u,d,A,c.beta,c.q,c.j);
                else
                    e=d-w.'*u;omega=1;
                    if method==2
                        z=c.gamma*e;
                        if abs(z)>=1e-12,omega=tanh(z)/z;end
                    end
                    r=u.'*P*u+c.r0/omega;
                    if ~isfinite(r) || r<=0,error('Invalid innovation denominator.');end
                    K=A*P*u/r;
                    w=A*w+K*e;
                    P=A*P*A.'-K*K.'*r+c.q*I;
                    P=(P+P.')/2;
                end
                if any(~isfinite(w)) || any(~isfinite(P(:)))
                    error('Nonfinite state or matrix.');
                end
                deviation(run,i)=(d-prediction)^2;
                if ~isfinite(deviation(run,i)),error('Nonfinite squared deviation.');end
            catch err
                failed(a)=failed(a)+1;
                if isempty(issue{a})
                    issue{a}=sprintf('run=%d sample=%d: %s',run,i,err.message);
                end
                break
            end
        end
    end
    if failed(a)==0
        curves(a,:)=10*log10(max(median(deviation,1),realmin));
        for win=1:3,metrics(a,win)=mean(curves(a,p.windows{win}));end
    end
end
end

function t=recovery(curve,post,change)
t=NaN;
if ~isfinite(post) || any(~isfinite(curve)),return;end
for i=change:numel(curve)-19
    if all(curve(i:i+19)<=post+3),t=i-change;return;end
end
end
