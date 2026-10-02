function bootstrapResults = bootstrap_breakpoint_comparison( ...
    dataByBiomarker, sigBiomarkers, B)

if nargin < 3
    B = 2000;
end
% bootstrap_breakpoint_comparison
% Implementation of bootstrap comparison of patient group-only breakpoint
% estimates between EEG theta/beta ratio and all other significant biomarkers.
%
% INPUT:
%   dataByBiomarker : struct with fields named after biomarkers.
%                     Each field contains a table with variables:
%                       - Biomarker (numeric, biomarker values)
%                       - Time (numeric, time relative to last
%                       follow-up/clinical onset)
%                       - Age (numeric)
%                       - Sex (categorical)
%                       - Group (categorical, where Group 1 is patient
%                       group)
%                       - SubjectID (categorical)
%   sigBiomarkers    : cell array with names of bioamrkers carried forward
%
% OUTPUT:
%   bootstrapResults : table containing
%       Biomarker_A, Biomarker_B
%       Tau_A, Tau_B
%       Delta_Tau
%       Delta_CI_low, Delta_CI_high
%       Tau_A_CI_low, Tau_A_CI_high
%       Tau_B_CI_low, Tau_B_CI_high
%       N_valid_bootstrap
%       Significant_difference

% -------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------
B = 2000;

idx_theta = find(strcmp(sigBiomarkers,'theta_beta'));
otherIdx  = setdiff(1:numel(sigBiomarkers), idx_theta);

pairs = [idx_theta * ones(numel(otherIdx),1), otherIdx(:)];
nPairs = size(pairs,1);

bootstrapResults = table();
rng(12345);

% -------------------------------------------------------------------------
% Loop through theta_beta vs each comparator biomarker
% -------------------------------------------------------------------------
for p = 1:nPairs

    biomarkerA = sigBiomarkers{pairs(p,1)};   % theta_beta
    biomarkerB = sigBiomarkers{pairs(p,2)};   % comparator

    dataA = dataByBiomarker.(biomarkerA);
    dataB = dataByBiomarker.(biomarkerB);

    % Patient subjects only
    dataA_conv = dataA(dataA.Group == 1,:);
    dataB_conv = dataB(dataB.Group == 1,:);

    subjA = unique(dataA_conv.SubjectID);
    subjB = unique(dataB_conv.SubjectID);

    commonSubjects = intersect(subjA,subjB);

    if numel(commonSubjects) < 5
        continue
    end

    nSubjects = numel(commonSubjects);

    % -------------------------------------------------------------------------
    % Preallocate bootstrap arrays
    % -------------------------------------------------------------------------
    tauA_boot = nan(B,1);
    tauB_boot = nan(B,1);
    deltaTau  = nan(B,1);
  
    bootError = strings(B,1);

    % ---------------------------------------------------------------------
    % Bootstrap (subject-level resampling)
    % ---------------------------------------------------------------------
    parfor b = 1:B

        sampledIdx = randi(nSubjects,nSubjects,1);
        sampledSubjects = commonSubjects(sampledIdx);

        bootA_rows = cell(nSubjects,1);
        bootB_rows = cell(nSubjects,1);

        for s = 1:nSubjects
            thisSubject = sampledSubjects(s);

            rowsA = dataA(dataA.SubjectID == thisSubject & dataA.Group == 1,:);
            rowsB = dataB(dataB.SubjectID == thisSubject & dataB.Group == 1,:);

            if isempty(rowsA) || isempty(rowsB)
                continue
            end

            newID = categorical({sprintf('boot_%d_%d',b,s)});
            rowsA.SubjectID = repmat(newID,height(rowsA),1);
            rowsB.SubjectID = repmat(newID,height(rowsB),1);

            bootA_rows{s} = rowsA;
            bootB_rows{s} = rowsB;
        end

        bootA = vertcat(bootA_rows{:});
        bootB = vertcat(bootB_rows{:});

        if isempty(bootA) || isempty(bootB)
            continue
        end

        try
            resultA = fit_piecewise_model_patients(bootA);
            resultB = fit_piecewise_model_patients(bootB);
        
            tauA_boot(b) = resultA.tau_hat;
            tauB_boot(b) = resultB.tau_hat;
            deltaTau(b)  = tauA_boot(b) - tauB_boot(b);
        
        catch ME
            bootError(b) = string(ME.message);
        end
    end

    % ---------------------------------------------------------------------
    % Clean bootstrap results
    % ---------------------------------------------------------------------
    valid = ~isnan(deltaTau);

    deltaValid = deltaTau(valid);
    tauA_valid = tauA_boot(valid);
    tauB_valid = tauB_boot(valid);

    nValid = numel(deltaValid);

    if nValid == 0
        continue
    end

    CI_delta = prctile(deltaValid,[2.5 97.5]);
    CI_tauA  = prctile(tauA_valid,[2.5 97.5]);
    CI_tauB  = prctile(tauB_valid,[2.5 97.5]);

    % ---------------------------------------------------------------------
    % Original (non-bootstrap) estimates
    % ---------------------------------------------------------------------
    originalA = fit_piecewise_model_patients(dataA);
    originalB = fit_piecewise_model_patients(dataB);

    tauA_original = originalA.tau_hat;
    tauB_original = originalB.tau_hat;

    delta_original = tauA_original - tauB_original;

    significantDifference = (CI_delta(1) > 0 || CI_delta(2) < 0);

    % ---------------------------------------------------------------------
    % Store results
    % ---------------------------------------------------------------------
    newRow = table( ...
        string(biomarkerA), ...
        string(biomarkerB), ...
        tauA_original, ...
        tauB_original, ...
        delta_original, ...
        CI_delta(1), CI_delta(2), ...
        CI_tauA(1), CI_tauA(2), ...
        CI_tauB(1), CI_tauB(2), ...
        nValid, ...
        significantDifference, ...
        'VariableNames',{ ...
            'Biomarker_A', ...
            'Biomarker_B', ...
            'Tau_A', ...
            'Tau_B', ...
            'Delta_Tau', ...
            'Delta_CI_low', ...
            'Delta_CI_high', ...
            'Tau_A_CI_low', ...
            'Tau_A_CI_high', ...
            'Tau_B_CI_low', ...
            'Tau_B_CI_high', ...
            'N_valid_bootstrap', ...
            'Significant_difference'});

    bootstrapResults = [bootstrapResults; newRow];

end

end
