function plot_piecewise_models(dataByBiomarker, biomarkerList)
% Plotting code for visualising z-scored (patinet-only)
% piecewise linear trajectories for multiple biomarkers.
%
% For each biomarker:
%   1. patient-only data are extracted.
%   2. Biomarker values are z-scored.
%   3. Breakpoint tau and post-break slope are estimated using
%      fit_piecewise_linear_converters.
%   4. A CI band and tau vertical line are plotted.

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
%   biomarkerList    ;  cell array of biomarker names to be plotted

% Colours
colorMap = lines(numel(biomarkerList));
colorStruct = containers.Map(biomarkerList, num2cell(colorMap,2));

% Time grid
t_range = linspace(-7,0,1000)';

% Storage
results = struct();
zslope  = struct();

% Fit models
for i = 1:numel(biomarkerList)

    biom = biomarkerList{i};
    D    = dataByBiomarker.(biom);

    D = D(D.Group == 1,:);
    if isempty(D)
        continue
    end

    mu = mean(D.Biomarker,'omitnan');
    sd = std(D.Biomarker,'omitnan');
    D.Z = (D.Biomarker - mu) ./ sd;

    Dfit = D;
    Dfit.Biomarker = D.Z;

    R = fit_piecewise_model_patients(Dfit);

    results.(biom) = R;
    zslope.(biom)  = R.post_slope;
end

% Plot
figure('Color','w'); hold on;
handles = gobjects(numel(biomarkerList),1);

for i = 1:numel(biomarkerList)

    biom = biomarkerList{i};
    if ~isfield(results, biom)
        continue
    end

    R = results.(biom);
    slope = abs(zslope.(biom));

    tau_hat = R.tau_hat;
    tau_CI  = R.tau_CI;

    c = colorStruct(biom);

    trajL = slope .* max(t_range - tau_CI(1), 0);
    trajU = slope .* max(t_range - tau_CI(2), 0);

    fill([t_range; flipud(t_range)], ...
         [trajL; flipud(trajU)], ...
         c, 'FaceAlpha',0.2, 'EdgeColor','none');

    handles(i) = xline(tau_hat,'--','Color',c,'LineWidth',2,'DisplayName',biom);
end

xlabel('Years to clinical onset');
ylabel('Z-scored biomarker value');
legend(handles,'Location','Best');
set(gca,'FontSize',14,'LineWidth',1.2,'TickDir','out','Box','off');

hold off;

end
