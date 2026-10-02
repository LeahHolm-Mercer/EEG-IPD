%%  EEG Spectral Identifiability Analysis 
%
%  This code:
%
%   1. Constructs feature vectors A and B for each EEG recording (20 epochs
%   in each)
%   2. Computes Pearson correlations between all A and B pairs → similarity matrix
%   3. Computes self‑rank for each EEG
%   4. Computes identity strength:
%        - z‑score correlations relative to between‑subject null
%   5. Mixed‑effects model:
%        Z ~ Condition + (1|ID_i) + (1|ID_j)
%   6. Longitudinal identity strength:
%        - Fisher z‑transform baseline→follow‑up correlations
%        - Linear slope per participant
%        - Group comparison via permutation test (10,000)
%        - Bootstrap CI for group mean slopes (5,000)
%   7. Plots:
%        - Condition distributions (Fig. 6A/B)
%        - Longitudinal trajectories (Fig. 6C)
%

clear; clc;
rng(12345);
% -------------------------------------------------------------------------
% 1. Load feature vectors and metadata
% -------------------------------------------------------------------------
%
% Feature vectors A and B contain the features derived from each 20 epochs from each recording
% Each EEG recording contributes 40 resting-state epochs.
% The epochs are randomly split into two sets of 20.
% Feature vectors A and B are generated from the two sets.
% Features are concatenated PSD from each scalp area
% 
% featureA   : matrix of feature vectors A (features x N recordings)
% featureB   : matrix of feature vectors B (features × N recordings)
% ageVec     : age per EEG recording (N x 1)
% sexVec     : sex per EEG recording (N x 1)
% groupVec   : diagnostic group per session (0 = Asymptomatic, 1 =
%              converter, 2 = patient), (N x 1)
% timeVec    : time relative to baseline recording for each individual per
%              recording (Nx1)
% IDVec      : subject ID per session (N x 1)

% -------------------------------------------------------------------------
% Diagnostic group coding
% -------------------------------------------------------------------------
% 0 = Control / asymptomatic
% 1 = Converter
% 2 = Symptomatic
%
% These numerical codes are retained for compatibility with the analysis
% pipeline. Meaningful categorical labels are created below.

GROUP_CONTROL     = 0;
GROUP_CONVERTER   = 1;
GROUP_SYMPTOMATIC = 2;

groupLabels = {'Control','Converter','Symptomatic'};

% Sort by diagnostic group (0,1,2) (for plotting)
desiredOrder = {'0','1','2'};
[~,sortIdx] = ismember(groupVec,desiredOrder);
[~,finalOrder] = sort(sortIdx);

featureA = featureA(:,finalOrder);
featureB = featureB(:,finalOrder);
ageVec   = ageVec(finalOrder);
sexVec   = sexVec(finalOrder);
groupVec = groupVec(finalOrder);
timeVec  = timeVec(finalOrder);
IDVec    = IDVec(finalOrder);

Group = categorical(groupVec, ...
    [GROUP_CONTROL GROUP_CONVERTER GROUP_SYMPTOMATIC], ...
    groupLabels);

Group = reordercats(Group,groupLabels);

% -------------------------------------------------------------------------
% 2. Compute Pearson correlation between all A and B recordign pairs
% -------------------------------------------------------------------------

featureA = zscore(featureA,0,2);
featureB = zscore(featureB,0,2);

nEEG = numel(groupVec);
correlations = nan(nEEG);

for i = 1:nEEG
    for j = 1:nEEG
        x = featureA(:,i);
        y = featureB(:,j);
        valid = ~isnan(x) & ~isnan(y);
        correlations(i,j) = corr(x(valid),y(valid));
    end
end

% -------------------------------------------------------------------------
% 3. Compute self‑rank for each EEG
% -------------------------------------------------------------------------

selfRank = nan(nEEG,1);
for i = 1:nEEG
    r_self = correlations(i,i);
    r_all  = correlations(i,:);
    [~,order] = sort(r_all,'descend');
    selfRank(i) = find(order==i);
end

% -------------------------------------------------------------------------
% 4. Compute identity strength: z‑score relative to between‑subject null
% -------------------------------------------------------------------------

upperMask = triu(true(nEEG),1);
betweenMask = (IDVec ~= IDVec') & upperMask;

betweenVals = correlations(betweenMask);
muB = mean(betweenVals);
sdB = std(betweenVals);

z_withinEEG = (diag(correlations) - muB) ./ sdB;

withinSubjMask = (IDVec == IDVec') & upperMask;
withinSubjVals = correlations(withinSubjMask);
z_withinSubj = (withinSubjVals - muB) ./ sdB;

z_between = (betweenVals - muB) ./ sdB;

% -------------------------------------------------------------------------
% 5. Mixed‑effects model: Z ~ Condition + (1|ID_i) + (1|ID_j)
% -------------------------------------------------------------------------

Condition = [
    repmat("WithinEEG",nEEG,1);
    repmat("Longitudinal",numel(z_withinSubj),1);
    repmat("Between",numel(z_between),1)
];

Condition = categorical(Condition, ...
    {'WithinEEG','Longitudinal','Between'});

Z_all = [z_withinEEG; z_withinSubj; z_between];

[iSubj,jSubj] = find(withinSubjMask);
[iBetween,jBetween] = find(betweenMask);

ID_i = [IDVec(:); IDVec(iSubj); IDVec(iBetween)];
ID_j = [IDVec(:); IDVec(jSubj); IDVec(jBetween)];

tblC = table(Z_all,categorical(Condition),categorical(ID_i),categorical(ID_j), ...
    'VariableNames',{'Z','Condition','ID_i','ID_j'});

lmeC = fitlme(tblC,'Z ~ Condition + (1|ID_i) + (1|ID_j)');
anovaC = anova(lmeC);
coefC = lmeC.Coefficients;

%% -------------------------------------------------------------------------
% Longitudinal identity strength (baseline → follow‑up)
% -------------------------------------------------------------------------

uniqueIDs = unique(IDVec);
EEGcol = [];
IDcol  = [];
DeltaTime = [];
Zcol = [];
Groupcol = [];

controlIdx = find(groupVec=='0');

for u = 1:numel(uniqueIDs)
    thisID = uniqueIDs(u);
    idx = find(IDVec==thisID);
    if numel(idx)<2, continue; end

    [~,order] = sort(timeVec(idx));
    idx = idx(order);

    baseIdx = idx(1);
    baseTime = timeVec(baseIdx);

    r_self = correlations(baseIdx,baseIdx);
    baseZ = atanh(r_self);

    EEGcol = [EEGcol; baseIdx];
    IDcol  = [IDcol; thisID];
    DeltaTime = [DeltaTime; 0];
    Zcol = [Zcol; baseZ];
    Groupcol = [Groupcol; groupVec(baseIdx)];

    for k = 2:numel(idx)
        laterIdx = idx(k);
        deltaT = timeVec(laterIdx) - baseTime;
        r_long = correlations(baseIdx,laterIdx);
        z_long = atanh(r_long);

        EEGcol = [EEGcol; laterIdx];
        IDcol  = [IDcol; thisID];
        DeltaTime = [DeltaTime; deltaT];
        Zcol = [Zcol; z_long];
        Groupcol = [Groupcol; groupVec(laterIdx)];
    end
end

tbl = table(EEGcol,categorical(IDcol),DeltaTime,Zcol,categorical(Groupcol), ...
    'VariableNames',{'EEG','ID','DeltaTime','Z','Group'});

% Linear slope per participant

uniqueIDs = unique(tbl.ID);
slope = nan(numel(uniqueIDs),1);
intercept = nan(numel(uniqueIDs),1);
groupSubj = categorical(nan(numel(uniqueIDs),1));

for i = 1:numel(uniqueIDs)
    thisID = uniqueIDs(i);
    rows = tbl(tbl.ID==thisID,:);
    if height(rows)<2, continue; end

    t = rows.DeltaTime;
    z = rows.Z;

    X = [ones(size(t)) t];
    b = X\z;

    intercept(i) = b(1);
    slope(i)     = b(2);
    groupSubj(i) = rows.Group(1);
end

valid = ~isnan(slope);
slope = slope(valid);
intercept = intercept(valid);
groupSubj = groupSubj(valid);
uniqueIDs = uniqueIDs(valid);

tblSlope = table(uniqueIDs,groupSubj,slope,intercept, ...
    'VariableNames',{'ID','Group','Slope','Intercept'});

% -------------------------------------------------------------------------
% Permutation tests for all pairwise group differences in slopes
% -------------------------------------------------------------------------

groupCodes = [GROUP_CONTROL GROUP_CONVERTER GROUP_SYMPTOMATIC];
groupNames = {'Control','Converter','Symptomatic'};

pairIdx = nchoosek(1:numel(groupCodes),2);
nPairs = size(pairIdx,1);

nPerm = 10000;

results = struct();
results.pairwisePermutation = struct();

pairwiseResults = table();

rng(12345);

for p = 1:nPairs

    gA = groupCodes(pairIdx(p,1));
    gB = groupCodes(pairIdx(p,2));

    nameA = groupNames{pairIdx(p,1)};
    nameB = groupNames{pairIdx(p,2)};

    sA = tblSlope.Slope(tblSlope.Group == gA);
    sB = tblSlope.Slope(tblSlope.Group == gB);

    sA = sA(~isnan(sA));
    sB = sB(~isnan(sB));

    if isempty(sA) || isempty(sB)
        warning('No valid slopes for %s vs %s.',nameA,nameB);
        continue;
    end

    observedDifference = mean(sA) - mean(sB);

    allSlopes = [sA; sB];

    groupLabelsPerm = [ ...
        zeros(numel(sA),1); ...
        ones(numel(sB),1)];

    permDifference = nan(nPerm,1);

    for k = 1:nPerm

        shuffledLabels = ...
            groupLabelsPerm(randperm(numel(groupLabelsPerm)));

        permDifference(k) = ...
            mean(allSlopes(shuffledLabels == 0)) ...
            - ...
            mean(allSlopes(shuffledLabels == 1));

    end

    pPermutation = ...
        (sum(abs(permDifference) >= abs(observedDifference)) + 1) ...
        / (nPerm + 1);

    newRow = table( ...
        string(nameA), ...
        string(nameB), ...
        numel(sA), ...
        numel(sB), ...
        observedDifference, ...
        pPermutation, ...
        nPerm, ...
        'VariableNames', ...
        {'Group1','Group2','N1','N2', ...
         'ObservedDifference','PermutationP','NPerm'});

    pairwiseResults = [pairwiseResults; newRow];

    results.pairwisePermutation(p).Group1 = label1;
    results.pairwisePermutation(p).Group2 = label2;
    results.pairwisePermutation(p).ObservedDifference = observedDifference;
    results.pairwisePermutation(p).P = pPermutation;
    results.pairwisePermutation(p).N_Group1 = n1;
    results.pairwisePermutation(p).N_Group2 = n2;
    results.pairwisePermutation(p).PermutationDistribution = permDifference;

end

disp(pairwiseResults);

% Bootstrap CI for group mean slopes (5,000) (reported)

nBoot = 5000;
meanSlope = zeros(numel(groups),1);
ciLower = zeros(numel(groups),1);
ciUpper = zeros(numel(groups),1);

for g = 1:numel(groups)
    s = tblSlope.Slope(tblSlope.Group==groups{g});
    meanSlope(g) = mean(s);
    bootstat = bootstrp(nBoot,@mean,s);
    ci = prctile(bootstat,[2.5 97.5]);
    ciLower(g) = ci(1);
    ciUpper(g) = ci(2);
end


% -------------------------------------------------------------------------
% 7. Plots
% -------------------------------------------------------------------------

% Plot similarity matrix (Fig. 6A)
nEEG = size(correlations,1);
figure('Color',[1 1 1]); 
imagesc(correlations);
cmap = [ ...
    0.18 0.31 0.55
    0.32 0.49 0.72
    0.55 0.67 0.82
    0.85 0.88 0.92
    0.96 0.85 0.80
    0.86 0.55 0.47
    0.70 0.26 0.22];
colormap(interp1(linspace(-1,1,size(cmap,1)), cmap, linspace(-1,1,256)));

caxis([-1 1])
cb = colorbar;
cb.Color = [0.2 0.2 0.2];

title('A–B Correlation Matrix','FontWeight','bold');
xlabel('Session B'); 
ylabel('Session A');

set(gca,'FontSize',12,...
        'LineWidth',1.2,...
        'Box','off',...
        'TickDir','out');

hold on;

groups = unique(groupVec,'stable');
startIdx = 1;

for g = 1:numel(groups)
    idx = find(groupVec == groups(g));
    nG  = numel(idx);
    
    rectangle('Position',[startIdx-0.5, startIdx-0.5, nG, nG], ...
        'EdgeColor',[0.3 0.3 0.3], ...   
        'LineWidth',1.8);
    
    startIdx = startIdx + nG;
end

hold off;

% Plot Fig. 6B

figure('Color','w'); hold on;
colors = [0.18 0.31 0.55; 0.75 0.35 0.28; 0.22 0.60 0.57];
conds = {"WithinEEG","Longitudinal","Between"};

for k = 1:3
    grp = conds{k};
    vals = Z_all(strcmp(Condition,grp));
    [f,xi] = ksdensity(vals);
    f = f./max(f);
    fill([xi fliplr(xi)], [f zeros(size(f))], colors(k,:), ...
        'FaceAlpha',0.35, 'EdgeColor',colors(k,:), 'LineWidth',1.2);
end

legend(conds,'Location','best');
xlabel('z-scored correlation');
ylabel('Normalised density');
title('Identity score distributions');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','off','TickDir','out');
hold off;


% Plot longitudinal trajectories (Fig. 6C)

figure('Color','w'); hold on;
colors = lines(numel(groups));

uniqueIDs = unique(tbl.ID);

for u = 1:numel(uniqueIDs)
    rows = tbl(tbl.ID==uniqueIDs(u),:);
    if height(rows)<2, continue; end
    g = find(strcmp(groups,string(rows.Group(1))));
    plot(rows.DeltaTime,rows.Z,'-o','Color',colors(g,:), ...
        'MarkerFaceColor',colors(g,:));
end

tgrid = linspace(min(tbl.DeltaTime),max(tbl.DeltaTime),100)';

for g = 1:numel(groups)
    m = meanSlope(g);
    lo = ciLower(g);
    hi = ciUpper(g);
    int = mean(tblSlope.Intercept(tblSlope.Group==groups{g}));

    y  = int + m*tgrid;
    yL = int + lo*tgrid;
    yU = int + hi*tgrid;

    fill([tgrid;flipud(tgrid)], [yL;flipud(yU)], colors(g,:), ...
        'FaceAlpha',0.15,'EdgeColor','none');
    plot(tgrid,y,'Color',colors(g,:),'LineWidth',3);
end

xlabel('Years from baseline');
ylabel('Fisher z-transformed similarity');
title('Longitudinal identity strength');
set(gca,'FontSize',14,'LineWidth',1.2,'Box','off','TickDir','out');
hold off;

% -------------------------------------------------------------------------
% 8. Store results
% -------------------------------------------------------------------------

results.similarity.correlations = correlations;
results.similarity.selfRank = selfRank;

results.identity.muBetween = muB;
results.identity.sdBetween = sdB;
results.identity.zWithinEEG = z_withinEEG;
results.identity.zWithinSubject = z_withinSubj;
results.identity.zBetween = z_between;

results.identity.LME = lmeC;
results.identity.LME_ANOVA = anova(lmeC);

results.longitudinal.table = tbl;
results.longitudinal.slopes = tblSlope;
results.longitudinal.observedDifference = obsDiff;
results.longitudinal.permutationDifference = permDiff;
results.longitudinal.pPermutation = p_perm;
results.longitudinal.nPermutation = nPerm;

results.longitudinal.meanSlope = meanSlope;
results.longitudinal.ciLower = ciLower;
results.longitudinal.ciUpper = ciUpper;
results.longitudinal.nBootstrap = nBoot;
