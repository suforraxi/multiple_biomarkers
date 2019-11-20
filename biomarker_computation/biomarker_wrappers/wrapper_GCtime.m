% wrapper GC time domain (similar to Park 2018)
function [ bio_vals, extra ] = wrapper_GCtime(cfg,data)

extra    = [];
bio_vals = [];



[~,fName,~]                  = fileparts(cfg.datasetName);
[ res_channel, artefact_T]   = get_metadata(cfg.inDir_data,fName);
[idxArtefact ,idx_art_trial] = find_artefacts_epochs(data.sampleinfo,data.label,artefact_T);

[~,outdata.hdr.datasetName,~] = fileparts(cfg.datasetName); 

cfgRedefine.trials = ~idx_art_trial;
data               = ft_redefinetrial(cfgRedefine,data);

ntrial = size(data.trial,2);

try
    res         = compute_GCtime(cfg,data);

    aux.m      = res.m;
    aux.pval   = res.pval;
    aux.sig    = res.sig;
    aux.morder = res.morder;

    extra{1}    = res;

    bio_vals{1} = (nansum(aux.m,1) / (size(aux.m,1)-1))'; % out strength
catch ME
    aux.m      = [];
    aux.pval   = [];
    aux.sig    = [];
    aux.morder = [];

    extra{1}    = [];

    bio_vals{1} = []; 

end




function res = compute_GCtime(cfg,a)
% a matrix (channels X time)


% detrend
%pp_a = detrend(a');



% first order diffrentiation

cfgDiff.differentiate = 'yes';

da = ft_preprocessing(cfgDiff,a);

[nChs,nSamples] = size(da.trial{1});
nTrials         = numel(da.trial);
 

data = zeros(nChs,nSamples,nTrials);

for i = 1 : numel(da.trial)
    
    data(:,:,i) = da.trial{i};
end

%zscore
z_a = zscore(data,[],2);


%% Parameters
ntrials   = nTrials;     % number of trials
nobs      = nSamples;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = cfg.momax;     % maximum model order for model order estimation

tstat     = '';     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')


X = z_a;

nvars = size(X,1);
%% Model order estimation

% Calculate information criteria up to max model order

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

[~,bmo_AIC] = min(AIC);
[~,bmo_BIC] = min(BIC);



% Plot information criteria.

% figure(1); clf;
% plot((1:momax)',[AIC BIC]);
% legend('AIC','BIC');

%amo = size(AT,3); % actual model order

fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
fprintf('best model order (BIC) = %d\n',bmo_BIC);
%fprintf('actual model order     = %d\n',amo);

% Select model order

if strcmpi(morder,'AIC')
    morder = bmo_AIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = bmo_BIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end

%% Granger causality estimation

% Calculate time-domain pairwise-conditional causalities. Return VAR parameters
% so we can check VAR.

ptic('\n*** GCCA_tsdata_to_pwcgc... ');
[F,A,SIG] = GCCA_tsdata_to_pwcgc(X,morder,regmode); % use same model order for reduced as for full regressions
ptoc;

% Check for failed (full) regression

assert(~isbad(A),'VAR estimation failed');

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Check VAR parameters (but don't bail out on error - GCCA mode is quite forgiving!)

rho = var_specrad(A);
fprintf('\nspectral radius = %f\n',rho);
if rho >= 1,       fprintf(2,'WARNING: unstable VAR (unit root)\n'); end
if ~isposdef(SIG), fprintf(2,'WARNING: residuals covariance matrix not positive-definite\n'); end

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat);
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.

% figure(2); clf;
% subplot(1,3,1);
% plot_pw(F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])
% 
% fprintf(2,'\nNOTE: no frequency-domain pairwise-conditional causality calculation in GCCA compatibility mode!\n');

res.m      = F;
res.pval   = pval;
res.sig    = sig;
res.morder = morder;