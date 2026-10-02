function result = fit_piecewise_model_patients(data)
% fit_piecewise_model_patients
% Implementation of breakpoint estimation for patient-only data.
%
% INPUT:
%   data : table with variables
%          - Biomarker (numeric, biomarker values)
%          - Time (numeric, relative to last follow-up/disease onset)
%          - Age (numeric)
%          - Sex (categorical)
%          - Group (categorical, where patients = 1)
%          - SubjectID (categorical)
%
% OUTPUT:
%   result : struct containing
%       .tau_hat (best breakpoint)
%       .tau_CI (95% confidence intervals around best breakpoint)
%       .pre_slope (slope before tau)
%       .post_slope (slope after tau)
%       .model

% -------------------------------------------------------------------------
% Restrict to patients
% -------------------------------------------------------------------------
varsRequired = {'Biomarker','Time','Age','Sex','Group','SubjectID'};
data = rmmissing(data(:,varsRequired));
data = data(data.Group == 1, :);
data.SubjectID = categorical(data.SubjectID);

% -------------------------------------------------------------------------
% tau grid
% -------------------------------------------------------------------------
tau_grid = linspace(quantile(data.Time,0.05), 0, 50);
logLik   = nan(length(tau_grid),1);

% -------------------------------------------------------------------------
% Profile likelihood
% -------------------------------------------------------------------------
for k = 1:length(tau_grid)
    tau = tau_grid(k);

    tmp = data;
    tmp.PostTime = max(tmp.Time - tau, 0);

    try
        if numel(unique(tmp.Sex)) == 1
            mdl_tmp = fitlme(tmp, ...
                'Biomarker ~ Age + Time + PostTime + (PostTime|SubjectID)', ...
                'FitMethod','ML', ...
                'DummyVarCoding','effects');
        else
            mdl_tmp = fitlme(tmp, ...
                'Biomarker ~ Age + Sex + Time + PostTime + (PostTime|SubjectID)', ...
                'FitMethod','ML', ...
                'DummyVarCoding','effects');
        end

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
% 95% CI
% -------------------------------------------------------------------------
LLmax = logLik(idx);
CI_idx = find(2*(LLmax - logLik) <= 3.84);

tau_CI = [tau_grid(min(CI_idx)), tau_grid(max(CI_idx))];

% -------------------------------------------------------------------------
% Final model at tau_hat
% -------------------------------------------------------------------------
data.PostTime = max(data.Time - tau_hat, 0);

if numel(unique(data.Sex)) == 1
    mdl = fitlme(data, ...
        'Biomarker ~ Age + Time + PostTime + (PostTime|SubjectID)', ...
        'FitMethod','REML', ...
        'DummyVarCoding','effects');
else
    mdl = fitlme(data, ...
        'Biomarker ~ Age + Sex + Time + PostTime + (PostTime|SubjectID)', ...
        'FitMethod','REML', ...
        'DummyVarCoding','effects');
end

coefTable = mdl.Coefficients;
names     = string(coefTable.Name);

pre_slope  = coefTable.Estimate(find(names == "Time"));
post_slope = coefTable.Estimate(find(names == "PostTime"));

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
result = struct();
result.tau_hat    = tau_hat;
result.tau_CI     = tau_CI;
result.pre_slope  = pre_slope;
result.post_slope = post_slope;
result.model      = mdl;

end
