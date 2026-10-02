% =========================================================================
% Requires MATLAB 2023a with Statistics and Machine Learning Toolbox
% =========================================================================

% =========================================================================
%  Baseline group differences (ANCOVA)
% =========================================================================
% Required inputs:
%   baseline : table with columns
%       Group (categorical)
%       Age   (numeric)
%       Sex   (categorical)
%       delta, theta, alpha, beta, gamma, theta_beta  (numeric EEG metrics)
%
% This section reproduces the baseline ANCOVA results reported in the paper.
% =========================================================================

vars = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'theta_beta'};
nVars  = numel(vars);
groups = categories(baseline.Patient);

p_ancova = nan(nVars,1);

F_omnibus       = nan(nVars,1);
df_group_all    = nan(nVars,1);
df_error_all    = nan(nVars,1);
eta_p2_omnibus  = nan(nVars,1);
eta_p2_low      = nan(nVars,1);
eta_p2_high     = nan(nVars,1);

for v = 1:nVars

    var = vars{v};

    y   = baseline.(var);
    g   = categorical(baseline.Patient);
    age = baseline.age;
    sex = categorical(baseline.Sex);
    Sens = categorical(baseline.FS);

    valid = ~isnan(y) & ~isnan(age) & ~isundefined(sex);

    tbl = table(y(valid), g(valid), age(valid), sex(valid),...
        'VariableNames', {'y','Group','Age','Sex'});

    mdl = fitlm(tbl, 'y ~ Group + Age + Sex');

    stats = anova(mdl,'components');

    row = strcmp(stats.Properties.RowNames,'Group');

    df_group = stats.DF(row);
    df_error = mdl.DFE;

    F_group = stats.F(row);
    p_group = stats.pValue(row);

    p_ancova(v) = p_group;

    eta_p2 = (F_group * df_group) / ...
             (F_group * df_group + df_error);

    [eta_L, eta_U] = partialEta2CI(F_group, df_group, df_error, 0.05);


    F_omnibus(v)      = F_group;
    df_group_all(v)   = df_group;
    df_error_all(v)   = df_error;
    eta_p2_omnibus(v) = eta_p2;
    eta_p2_low(v)     = eta_L;
    eta_p2_high(v)    = eta_U;

    R2     = mdl.Rsquared.Ordinary;
    R2_adj = mdl.Rsquared.Adjusted;

    coefNames = mdl.Coefficients.Properties.RowNames;

    group1Row = contains(coefNames,'Group_1');
    group2Row = contains(coefNames,'Group_2');

    if sum(group1Row) ~= 1 || sum(group2Row) ~= 1
        error('Could not uniquely identify Group_1 and Group_2 coefficients.');
    end


    beta1 = mdl.Coefficients.Estimate(group1Row);
    SE1   = mdl.Coefficients.SE(group1Row);

    beta2 = mdl.Coefficients.Estimate(group2Row);
    SE2   = mdl.Coefficients.SE(group2Row);

    t1 = mdl.Coefficients.tStat(group1Row);
    t2 = mdl.Coefficients.tStat(group2Row);

    tcrit = tinv(0.975, df_error);

    CI95_beta1 = beta1 + [-tcrit, tcrit] * SE1;
    CI95_beta2 = beta2 + [-tcrit, tcrit] * SE2;

    eta_p2_group1 = (t1^2) / (t1^2 + df_error);
    eta_p2_group2 = (t2^2) / (t2^2 + df_error);

    [eta1_L, eta1_U] = partialEta2CI(t1^2, 1, df_error, 0.05);
    [eta2_L, eta2_U] = partialEta2CI(t2^2, 1, df_error, 0.05);

  
    fprintf('\n%s:\n', var);

    fprintf(['  Omnibus Group: F(%d,%d)=%.2f, p=%.4g, ', ...
             'partial eta^2=%.3f, 95%% CI [%.3f, %.3f], ', ...
             'R^2=%.3f, R^2adj=%.3f\n'], ...
        df_group, df_error, F_group, p_group, ...
        eta_p2, eta_L, eta_U, R2, R2_adj);

    fprintf(['  Group 1 vs 0: beta=%.3f, 95%% CI [%.3f, %.3f], ', ...
             'partial eta^2=%.3f, 95%% CI [%.3f, %.3f]\n'], ...
        beta1, CI95_beta1(1), CI95_beta1(2), ...
        eta_p2_group1, eta1_L, eta1_U);

    fprintf(['  Group 2 vs 0: beta=%.3f, 95%% CI [%.3f, %.3f], ', ...
             'partial eta^2=%.3f, 95%% CI [%.3f, %.3f]\n'], ...
        beta2, CI95_beta2(1), CI95_beta2(2), ...
        eta_p2_group2, eta2_L, eta2_U);

end


% FDR correction 

p_ancova_fdr = bh_fdr(p_ancova);

fprintf('\n=== FDR-adjusted ANCOVA p-values ===\n');

for v = 1:nVars
    fprintf('%s: q = %.4g\n', vars{v}, p_ancova_fdr(v));
end


% Pairwise comparisons

pairs  = nchoosek(1:numel(groups), 2);
nPairs = size(pairs,1);

pairwise_p     = nan(nVars, nPairs);
pairwise_p_fdr = nan(nVars, nPairs);

pairwise_eta_p2      = nan(nVars, nPairs);
pairwise_eta_p2_low  = nan(nVars, nPairs);
pairwise_eta_p2_high = nan(nVars, nPairs);

pairwise_beta    = nan(nVars, nPairs);
pairwise_beta_low  = nan(nVars, nPairs);
pairwise_beta_high = nan(nVars, nPairs);

for v = 1:nVars

    if p_ancova_fdr(v) > 0.05
        continue;
    end

    y   = baseline.(vars{v});
    g   = categorical(baseline.Patient);
    age = baseline.age;
    sex = categorical(baseline.Sex);
    Sens = categorical(baseline.MentalActivation);
    
    valid = ~isnan(y) & ~isnan(age) & ...
            ~isundefined(sex) & ~isundefined(Sens);
    
    tbl = table(y(valid), g(valid), age(valid), sex(valid), Sens(valid), ...
        'VariableNames', {'y','Group','Age','Sex','Sens'});
    
    lm = fitlm(tbl, 'y ~ Group + Age + Sex + Sens');

    coefNames = lm.CoefficientNames;
    df_error  = lm.DFE;

    for p = 1:nPairs

        g1 = groups{pairs(p,1)};
        g2 = groups{pairs(p,2)};


        c = zeros(1, numel(coefNames));

        if ~strcmp(g1, groups{1})
            c(strcmp(coefNames, ['Group_' g1])) = 1;
        end

        if ~strcmp(g2, groups{1})
            c(strcmp(coefNames, ['Group_' g2])) = -1;
        end

        [p_pair, F_pair] = coefTest(lm, c);

        pairwise_p(v,p) = p_pair;

        eta_pair = F_pair / (F_pair + df_error);

        pairwise_eta_p2(v,p) = eta_pair;

        [eta_L, eta_U] = partialEta2CI(F_pair, 1, df_error, 0.05);
        
        pairwise_eta_p2_low(v,p)  = eta_L;
        pairwise_eta_p2_high(v,p) = eta_U;

        beta_pair = c * lm.Coefficients.Estimate;

        se_pair = sqrt(c * lm.CoefficientCovariance * c');

        t_pair = beta_pair / se_pair;

        tcrit = tinv(0.975, df_error);

        pairwise_beta(v,p) = beta_pair;
        pairwise_beta_low(v,p) = beta_pair - tcrit * se_pair;
        pairwise_beta_high(v,p) = beta_pair + tcrit * se_pair;

    end

    pairwise_p_fdr(v,:) = bh_fdr(pairwise_p(v,:)')';

end

fprintf('\n=== Pairwise Comparisons ===\n');

for v = 1:nVars

    if p_ancova_fdr(v) > 0.05
        continue;
    end

    fprintf('\nVariable: %s\n', vars{v});

    fprintf(['Group1\tGroup2\tp_adj\t', ...
             'partial eta^2\t95%% CI eta^2\n']);

    for p = 1:nPairs

        if isnan(pairwise_p_fdr(v,p))
            continue;
        end

        g1 = groups{pairs(p,1)};
        g2 = groups{pairs(p,2)};

        fprintf('%s\t%s\t%.4g\t%.3f\t[%.3f, %.3f]\n', ...
            g1, g2, ...
            pairwise_p_fdr(v,p), ...
            pairwise_eta_p2(v,p), ...
            pairwise_eta_p2_low(v,p), ...
            pairwise_eta_p2_high(v,p));
    end
end

% Violin plots
nGroups = numel(groups);
grpColors = lines(nGroups);   

for v = 1:numel(vars)

    figure; hold on;

    x = baseline.(vars{v});
    width = 0.35;
   
    for iG = 1:nGroups
        grpData = x(g == groups{iG});

        simple_violin(grpData, iG, width, grpColors(iG,:));

        jitter = (rand(size(grpData)) - 0.5) * 0.15;
        scatter(iG + jitter, grpData, ...
            20, grpColors(iG,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.7, ...
            'MarkerEdgeColor', 'none');
    end
    set(gca, 'XTick', [], 'XTickLabel', []);
end


% =========================================================================
% Qualitative EEG analysis: logistic regression
% =========================================================================
% EEG clinical reports were classified as:
%   0 = within normal limits
%   1 = abnormal or non-specific irregularities of uncertain significance
%
% The binary EEG classification was modelled using logistic regression,
% with patient group as the predictor. The asymptomatic group was used
% as the reference category. Group-specific contrasts were assessed using
% Wald tests.
% Requires 'baseline' table as for the ANCOVA.
% =========================================================================

EEGtxt = baseline.EEGComments;

% NOTE:
% This assumes that a non-empty EEG comment corresponds to an abnormal
% or non-specific EEG report. 

EEG_abnormal = ~cellfun(@isempty, EEGtxt);

tblQual = table( ...
    categorical(baseline.Group), ...
    double(EEG_abnormal), ...
    'VariableNames', {'Group','EEG_abnormal'});


groupsQual = categories(tblQual.Group);


tblQual.Group = reordercats(tblQual.Group, ...
    [{'Asymptomatic'}, groupsQual(~strcmp(groupsQual,'Asymptomatic'))]);

mdl_log = fitglm(tblQual, ...
    'EEG_abnormal ~ Group', ...
    'Distribution','binomial', ...
    'Link','logit');

coef = mdl_log.Coefficients;

disp(' ');
disp('Model coefficients and Wald tests:');
disp(coef);

coefNames = mdl_log.CoefficientNames;

groupCoefIdx = startsWith(coefNames,'Group_');

C_group = zeros(sum(groupCoefIdx), numel(coefNames));
C_group(:,groupCoefIdx) = eye(sum(groupCoefIdx));

p_group_overall = coefTest(mdl_log, C_group);

nonRefGroups = groupsQual(~strcmp(groupsQual,'Asymptomatic'));

p_contrast = NaN(numel(nonRefGroups),1);
beta_contrast = NaN(numel(nonRefGroups),1);
SE_contrast = NaN(numel(nonRefGroups),1);
OR_contrast = NaN(numel(nonRefGroups),1);
CI_low = NaN(numel(nonRefGroups),1);
CI_high = NaN(numel(nonRefGroups),1);

for i = 1:numel(nonRefGroups)

    groupName = nonRefGroups{i};

    coefIdx = strcmp(coefNames, ['Group_' groupName]);

 
    C = zeros(1,numel(coefNames));
    C(coefIdx) = 1;

    p_contrast(i) = coefTest(mdl_log, C);

    beta_contrast(i) = coef.Estimate(coefIdx);
    SE_contrast(i) = coef.SE(coefIdx);

    OR_contrast(i) = exp(beta_contrast(i));
    CI_low(i) = exp(beta_contrast(i) - 1.96*SE_contrast(i));
    CI_high(i) = exp(beta_contrast(i) + 1.96*SE_contrast(i));

    fprintf('\n%s vs Asymptomatic:\n', groupName);
    fprintf('  Beta = %.4f\n', beta_contrast(i));
    fprintf('  SE   = %.4f\n', SE_contrast(i));
    fprintf('  OR   = %.3f\n', OR_contrast(i));
    fprintf('  95%% CI = [%.3f, %.3f]\n', CI_low(i), CI_high(i));
    fprintf('  Wald p = %.6g\n', p_contrast(i));

end



% =========================================================================
%  Longitudinal change in converters (LME)
% =========================================================================
% Required input:
%   cData : table with columns
%       SubjectID (categorical)
%       Time      (numeric, years relative to onset)
%       Age       (numeric)
%       Sex       (categorical)
%       EEG metrics: delta, theta, alpha, beta, gamma, theta_beta
%
% This section reproduces the longitudinal LME results reported in the paper.
% =========================================================================

timeVar  = cData.time;
subjects = categorical(cData.SubjectID);
age      = cData.age;
sex      = categorical(cData.Sex);

nVars = numel(vars);

all_p = nan(nVars,1);
all_stats = cell(nVars,1);

for v = 1:nVars

    y = cData.(vars{v});

    tbl = table(y, timeVar, subjects, age, sex, ...
        'VariableNames', {'y','Time','Subject','Age','Sex'});

    lme = fitlme(tbl, 'y ~ Time + Age + Sex + (1|Subject)');

    coef = lme.Coefficients;
    row  = strcmp(coef.Name,'Time');

    beta = coef.Estimate(row);
    SE   = coef.SE(row);
    tval = coef.tStat(row);
    pval = coef.pValue(row);
    CI   = coef.Lower(row) + [0, coef.Upper(row)-coef.Lower(row)];

    all_p(v) = pval;

    all_stats{v} = struct( ...
        'Variable', vars{v}, ...
        'Beta', beta, ...
        'SE', SE, ...
        't', tval, ...
        'p', pval, ...
        'CI', CI ...
    );

    % Plot
    figure; hold on;

    uSubs = unique(subjects);
    for s = 1:numel(uSubs)
        idx = subjects == uSubs(s);
        plot(timeVar(idx), y(idx), '-o', ...
            'Color', [0 0 0 0.3], ...
            'LineWidth', 1.5, ...              
            'MarkerSize', 4, ...             
            'MarkerFaceColor', 'none', ... 
            'MarkerEdgeColor', 'k');
    end

    tGrid = linspace(min(timeVar), max(timeVar), 200)';
    
    subjDummy = repmat(subjects(1), length(tGrid), 1);
    ageMean   = mean(age);
    sexMode   = mode(sex);
    yMean     = mean(timeVar);
    
    newTbl = table( ...
        yMean * ones(length(tGrid),1), ...   % y
        tGrid, ...                           % Time
        subjDummy, ...                       % Subject
        repmat(ageMean,length(tGrid),1), ... % Age
        repmat(sexMode,length(tGrid),1), ... % Sex
        'VariableNames', {'y','Time','Subject','Age','Sex'} ...
    );
    
    [yPred, yCI] = predict(lme, newTbl, ...
        'Conditional', false, ...
        'Alpha', 0.05);

    col = colors(2,:);
    patch([tGrid; flipud(tGrid)], ...
          [yCI(:,1); flipud(yCI(:,2))], ...
          col, 'FaceAlpha',0.2, 'EdgeColor','none');

    plot(tGrid, yPred, 'Color', col, 'LineWidth', 2);

    grid on;
    hold off;

end

adj_p = bh_fdr(all_p);

% =========================================================================
%  Associations between EEG metrics and other biomarkers (LME)
% =========================================================================
% Required input:
%   cData must contain columns for:
%       NfL, GFAP, PIQ, Stroop, and other biomarkers wanted to compare
%
% This section reproduces the EEG–biomarker association results reported.
% =========================================================================

vars = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'theta_beta'};
Compvars = {'NfL', 'GFAP', 'PIQ', 'Stroop', 'TMT-A', 'TMT-B'};
nComp = numel(Compvars);
nEEG = numel(vars);
colors = lines(nGroups);

% Store results
beta_mat = nan(nEEG, nComp);
se_mat   = nan(nEEG, nComp);
t_mat    = nan(nEEG, nComp);
p_mat    = nan(nEEG, nComp);
ci_low   = nan(nEEG, nComp);
ci_high  = nan(nEEG, nComp);

for v = 1:nEEG
    for n = 1:nComp

        x = cData.(Compvars{n});
        nan_idx = isnan(x);
        
        currData = cData(~nan_idx, :);
        x = x(~nan_idx);
        y = currData.(vars{v});
       
        tbl = table(y, x, currData.age, categorical(currData.CohortID), ...
                    categorical(currData.Sex), ...
            'VariableNames', {'y','x','Age','Subject','Sex'});

        % Fit model: EEG ~ NP + Age + Sex + (1|Subject)
        lme = fitlme(tbl, 'y ~ x + Age + Sex + (1|Subject)');

        % Extract NP coefficient
        coef = lme.Coefficients;
        row = strcmp(coef.Name, 'x');

        beta_mat(v,n)  = coef.Estimate(row);
        se_mat(v,n)    = coef.SE(row);
        t_mat(v,n)     = coef.tStat(row);
        p_mat(v,n)     = coef.pValue(row);
        ci_low(v,n)    = coef.Lower(row);
        ci_high(v,n)   = coef.Upper(row);

        % Plot
        figure; hold on;

        scatter(tbl.x, tbl.y, 25, 'k', 'filled', 'MarkerFaceAlpha', 0.4);

        % Prediction grid
        gGrid = linspace(min(tbl.x), max(tbl.x), 100)';
        ageMean = mean(tbl.Age, 'omitnan');
        sexMode = mode(tbl.Sex);
        subjDummy = repmat(tbl.Subject(1), length(gGrid), 1);

        newTbl = table(gGrid, repmat(ageMean,length(gGrid),1), subjDummy, ...
                       repmat(sexMode,length(gGrid),1), ...
                       'VariableNames', {'x','Age','Subject','Sex'});

        [yPred, yCI] = predict(lme, newTbl, ...
            'Conditional', false, ...
            'Alpha', 0.05);

        patch([gGrid; flipud(gGrid)], ...
              [yCI(:,1); flipud(yCI(:,2))], ...
              [0.8 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.3);

        plot(gGrid, yPred, 'b-', 'LineWidth', 2);

        hold off;
    end
end

% FDR correction across all comparisons
p_vec = p_mat(:);
p_fdr_vec = bh_fdr(p_vec);
p_fdr = reshape(p_fdr_vec, size(p_mat));
