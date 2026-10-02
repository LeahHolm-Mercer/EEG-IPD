%% demo_piecewise_analysis
% Synthetic demonstration of the piecewise-linear analysis pipeline.
%
% This script:
%   1. Generates synthetic longitudinal patient/control data.
%   2. Stores the data in dataByBiomarker with the same field names used
%      by the analysis functions.
%   3. Fits patient + control piecewise models.
%   4. Fits patient-only models for the significant biomarkers.
%   5. Runs the patient-only bootstrap breakpoint comparison.
%   6. Produces a plot of the synthetic data and fitted trajectories.
%
% Group coding:
%   0 = control
%   1 = patient/converter
%
% Time:
%   Years relative to disease onset, with onset at Time = 0.
%   Negative values therefore represent years before onset.

clear;
clc;
close all;

%% ------------------------------------------------------------------------
%  SETTINGS
% -------------------------------------------------------------------------

rng(12345);

% Number of subjects
nControls = 20;
nPatients = 20;

% Repeated measurements per subject
nTimePoints = 8;

% Synthetic biomarker names.
% These names must match the fields expected by the bootstrap function.
sigBiomarkers = {'NfL','GFAP','PIQ','Stroop','TMTB','theta_beta'};

% True synthetic breakpoints.
% These are deliberately different so that the breakpoint estimation
% and bootstrap comparison have something meaningful to detect.
trueTau.NfL        = -3.5;
trueTau.GFAP       = -2.5;
trueTau.PIQ        = -4.0;
trueTau.Stroop     = -3.0;
trueTau.TMTB       = -2.5;
trueTau.theta_beta = -2.0;

% Patient pre-breakpoint slopes
preSlope.NfL        = 0.02;
preSlope.GFAP       = 0.01;
preSlope.PIQ        = -0.02;
preSlope.Stroop     = 0.03;
preSlope.TMTB       = 0.04;
preSlope.theta_beta = 0.01;

% Additional slope after the breakpoint.
postSlope.NfL        = 0.25;
postSlope.GFAP       = 0.18;
postSlope.PIQ        = -0.20;
postSlope.Stroop     = 0.30;
postSlope.TMTB       = 0.35;
postSlope.theta_beta = 0.28;

% Measurement noise
noiseSD.NfL        = 0.20;
noiseSD.GFAP       = 0.18;
noiseSD.PIQ        = 0.20;
noiseSD.Stroop     = 0.25;
noiseSD.TMTB       = 0.25;
noiseSD.theta_beta = 0.15;

%% ------------------------------------------------------------------------
%  GENERATE SUBJECT INFORMATION
% -------------------------------------------------------------------------

nSubjects = nControls + nPatients;

SubjectID = categorical( ...
    arrayfun(@(x) sprintf('S%03d',x), ...
    1:nSubjects, 'UniformOutput', false)');

Group = [zeros(nControls,1); ones(nPatients,1)];

% Give each subject a fixed age and sex.
Age = 35 + 15*rand(nSubjects,1);

Sex = categorical( ...
    randi([0 1],nSubjects,1), ...
    [0 1], ...
    {'Female','Male'});

%% ------------------------------------------------------------------------
%  GENERATE LONGITUDINAL DATA FOR EACH BIOMARKER
% -------------------------------------------------------------------------

dataByBiomarker = struct();

for v = 1:numel(sigBiomarkers)

    biomarker = sigBiomarkers{v};

    rows = cell(nSubjects,1);

    for s = 1:nSubjects

        % Slightly jitter measurement times between subjects.
        time = linspace(-7,0,nTimePoints)' ...
             + 0.08*randn(nTimePoints,1);

        % Keep all observations within the analysis window.
        time = max(min(time,0),-7);
        time = sort(time);

        thisAge = Age(s);
        thisSex = repmat(Sex(s),nTimePoints,1);
        thisGroup = repmat(Group(s),nTimePoints,1);
        thisSubject = repmat(SubjectID(s),nTimePoints,1);

        tau = trueTau.(biomarker);

                % Baseline level (synthetic)
        switch biomarker
            case 'NfL'
                baseline = 1.0;
            case 'GFAP'
                baseline = 0.8;
            case 'PIQ'
                baseline = 100;
            case 'Stroop'
                baseline = 40;
            case 'TMTB'
                baseline = 70;
            case 'theta_beta'
                baseline = 0.5;
        end

        % Age effect (small)
        switch biomarker
            case 'NfL'
                ageEff = 0.005;
            case 'GFAP'
                ageEff = 0.004;
            case 'PIQ'
                ageEff = -0.10;
            case 'Stroop'
                ageEff = 0.06;
            case 'TMTB'
                ageEff = 0.08;
            case 'theta_beta'
                ageEff = 0.003;
        end
        ageEffect = ageEff * (thisAge - mean(Age));

        % Sex effect (small)
        switch biomarker
            case 'NfL'
                sexEff = 0.05;
            case 'GFAP'
                sexEff = 0.04;
            case 'PIQ'
                sexEff = 2;
            case 'Stroop'
                sexEff = 1.5;
            case 'TMTB'
                sexEff = 2.0;
            case 'theta_beta'
                sexEff = 0.02;
        end
        sexCode = double(thisSex == 'Male');
        sexEffect = sexEff * sexCode;

        % Control background slope
        switch biomarker
            case 'NfL'
                controlSlope = 0.01;
            case 'GFAP'
                controlSlope = 0.008;
            case 'PIQ'
                controlSlope = -0.01;
            case 'Stroop'
                controlSlope = 0.02;
            case 'TMTB'
                controlSlope = 0.03;
            case 'theta_beta'
                controlSlope = 0.005;
        end

        % Piecewise trajectory:
        %
        % pre-breakpoint:
        %     baseline + preSlope * Time
        %
        % post-breakpoint:
        %     baseline + preSlope * Time ...
        %       + postSlope * PostTime
        %
        % where PostTime = max(Time - tau, 0).
        postTime = max(time - tau,0);

        patientEffect = preSlope.(biomarker) .* time ...
                      + postSlope.(biomarker) .* postTime;

        controlTrajectory = controlSlope .* time;

        % Use the disease-associated piecewise trajectory only for patients.
        trajectory = zeros(nTimePoints,1);

        patientMask = (thisGroup == 1);
        controlMask = (thisGroup == 0);

        trajectory(patientMask) = patientEffect(patientMask);
        trajectory(controlMask) = controlTrajectory(controlMask);

        % Add measurement noise.
        y = baseline ...
          + ageEffect ...
          + sexEffect ...
          + trajectory ...
          + noiseSD.(biomarker).*randn(nTimePoints,1);

        rows{s} = table( ...
            thisSubject, ...
            thisGroup, ...
            repmat(thisAge,nTimePoints,1), ...
            thisSex, ...
            time, ...
            y, ...
            'VariableNames', ...
            {'SubjectID','Group','Age','Sex','Time','Biomarker'});
    end

    dataByBiomarker.(biomarker) = vertcat(rows{:});
end

%% ------------------------------------------------------------------------
%  CHECK THE GENERATED DATA
% -------------------------------------------------------------------------

disp('Synthetic data generated.');

for v = 1:numel(sigBiomarkers)

    biomarker = sigBiomarkers{v};
    D = dataByBiomarker.(biomarker);

    fprintf('%s: %d rows, %d subjects\n', ...
        biomarker, height(D), numel(unique(D.SubjectID)));
end

%% ------------------------------------------------------------------------
%  PATIENT + CONTROL PIECEWISE MODELS
% -------------------------------------------------------------------------

fprintf('\n=============================================\n');
fprintf('PATIENT + CONTROL PIECEWISE MODELS\n');
fprintf('=============================================\n');

modelResults = struct();

for v = 1:numel(sigBiomarkers)

    biomarker = sigBiomarkers{v};

    fprintf('\nFitting %s...\n', biomarker);

    D = dataByBiomarker.(biomarker);

    result = fit_piecewise_model(D);

    modelResults.(biomarker) = result;

    fprintf('  Estimated breakpoint: %.3f years\n', ...
        result.tau_hat);

    fprintf('  95%% CI: [%.3f, %.3f]\n', ...
        result.tau_CI(1), result.tau_CI(2));

    fprintf('  Control post-break slope: %.4f\n', ...
        result.group0_slope);

    fprintf('  Patient post-break slope: %.4f\n', ...
        result.group1_slope);

    fprintf('  Patient-control slope difference: %.4f\n', ...
        result.slope_difference);

    fprintf('  Interaction p-value: %.4g\n', ...
        result.p_two_sided);
end

%% ------------------------------------------------------------------------
%  FDR on interaction p-values → create sigBiomarkers
% -------------------------------------------------------------------------

allBiomarkers = fieldnames(modelResults);

allP = nan(numel(allBiomarkers),1);
for i = 1:numel(allBiomarkers)
    bm = allBiomarkers{i};
    allP(i) = modelResults.(bm).p_two_sided;
end

qvals = bh_fdr(allP);

sigMask = qvals < 0.05;
sigBiomarkers = allBiomarkers(sigMask);

fprintf('\nSignificant biomarkers after FDR:\n');
disp(sigBiomarkers);


%% ------------------------------------------------------------------------
%  PATIENT-ONLY MODELS
% -------------------------------------------------------------------------

fprintf('\n=============================================\n');
fprintf('PATIENT-ONLY BREAKPOINT MODELS\n');
fprintf('=============================================\n');

patientResults = struct();

for v = 1:numel(sigBiomarkers)

    biomarker = sigBiomarkers{v};

    fprintf('\nFitting patient-only model: %s...\n', biomarker);

    D = dataByBiomarker.(biomarker);

    result = fit_piecewise_model_patients(D);

    patientResults.(biomarker) = result;

    fprintf('  Estimated breakpoint: %.3f years\n', ...
        result.tau_hat);

    fprintf('  95%% CI: [%.3f, %.3f]\n', ...
        result.tau_CI(1), result.tau_CI(2));

    fprintf('  Pre-break slope: %.4f\n', ...
        result.pre_slope);

    fprintf('  Post-break slope: %.4f\n', ...
        result.post_slope);
end

%% ------------------------------------------------------------------------
%  BOOTSTRAP COMPARISON OF BREAKPOINTS
% -------------------------------------------------------------------------
%
% For a quick demonstration, 200 bootstrap samples (B) are recommended.
% For the final analysis, use B = 2000.

fprintf('\n=============================================\n');
fprintf('BOOTSTRAP BREAKPOINT COMPARISON\n');
fprintf('=============================================\n');

B = 2000;

try

    bootstrapResults = ...
        bootstrap_breakpoint_comparison( ...
            dataByBiomarker, sigBiomarkers, B);

    disp(bootstrapResults);

catch ME

    fprintf('\nBootstrap could not be run.\n');
    fprintf('Error message:\n%s\n\n',ME.message);

    fprintf(['Check that your bootstrap function has the updated ', ...
             'input arguments:\n']);
    fprintf('bootstrap_breakpoint_comparison(dataByBiomarker, sigBiomarkers, B)\n');

    bootstrapResults = table();
end

%% ------------------------------------------------------------------------
%  PLOT BREAKPOINTS
% -------------------------------------------------------------------------

plot_piecewise_models(dataByBiomarker, sigBiomarkers);
