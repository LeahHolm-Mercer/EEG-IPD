function result = fit_piecewise_model(data)
% fit_piecewise_model
% Implementation of breakpoint estimation using
% profile likelihood and mixed-effects modelling.
%
% INPUT:
%   data : table with variables
%          - Biomarker (numeric, value)
%          - Time (numeric, relative to last follow up or conversion)
%          - Age (numeric)
%          - Sex (categorical or char)
%          - Group (categorical, patient = 1, control = 0)
%          - SubjectID (categorical)
%
% OUTPUT:
%   result : struct containing
%       .tau_hat (breakpoint)
%       .tau_CI (95% CI)
%       .group0_slope (controls)
%       .group1_slope (patients)
%       .slope_difference (between controls and patients after tau_hat)
%       .p_two_sided
%       .model

% -------------------------------------------------------------------------
% Prepare variables
% -------------------------------------------------------------------------
data.SubjectID = categorical(data.SubjectID);
data.Group = categorical(data.Group);
data.Group = reordercats(data.Group, {'0','1'}); % make sure that control is the reference
varsRequired = {'Biomarker','Time','Age','Sex','Group','SubjectID'};
data = rmmissing(data(:,varsRequired));
% -------------------------------------------------------------------------
% Define tau search grid
% -------------------------------------------------------------------------
tau_grid = linspace(quantile(data.Time,0.05), 0, 50);
logLik   = nan(length(tau_grid),1);

% -------------------------------------------------------------------------
% Profile likelihood loop
% -------------------------------------------------------------------------
for k = 1:length(tau_grid)
    tau = tau_grid(k);

    tmp = data;
    tmp.PostTime = max(tmp.Time - tau, 0);

    try
        mdl_tmp = fitlme(tmp, ...
            'Biomarker ~ Age + Sex + PostTime + Group:PostTime + (PostTime|SubjectID)', ...
            'FitMethod','ML');

        logLik(k) = mdl_tmp.LogLikelihood;

    catch
        logLik(k) = NaN;
    end
end

% -------------------------------------------------------------------------
% Best tau
% -------------------------------------------------------------------------
[~, idx] = max(logLik);
tau_hat = tau_grid(idx);


% -------------------------------------------------------------------------
% 95% CI using profile likelihood
% -------------------------------------------------------------------------
LLmax = logLik(idx);
cutoff = LLmax - 0.5 * chi2inv(0.95,1);

inside = logLik >= cutoff;

if any(inside)
    tau_CI = [min(tau_grid(inside)), max(tau_grid(inside))];
else
    tau_CI = [NaN NaN];
end

% -------------------------------------------------------------------------
% Final model at tau_hat
% -------------------------------------------------------------------------
data.PostTime = max(data.Time - tau_hat, 0);

mdl = fitlme(data, ...
    'Biomarker ~ Age + Sex + PostTime + Group:PostTime + (PostTime|SubjectID)', ...
    'FitMethod','REML');

coefTable = mdl.Coefficients;
names     = string(coefTable.Name);

% -------------------------------------------------------------------------
% Extract slopes
% -------------------------------------------------------------------------
idx_post  = find(names == "PostTime");
beta0     = coefTable.Estimate(idx_post);

idx_group = find(contains(names,"Group") & contains(names,"PostTime"));
beta_grp  = coefTable.Estimate(idx_group);

group0_slope = beta0;
group1_slope = beta0 + beta_grp;

% -------------------------------------------------------------------------
% p-values
% -------------------------------------------------------------------------
p_two_sided = coefTable.pValue(idx_group);

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
result = struct();
result.tau_hat         = tau_hat;
result.tau_CI          = tau_CI;
result.group0_slope    = group0_slope;
result.group1_slope    = group1_slope;
result.slope_difference = beta_grp;
result.p_two_sided     = p_two_sided;
result.model           = mdl;

end
