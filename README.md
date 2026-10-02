# EEG Identifiability, Piecewise Modelling, and Statistical Analysis

This repository contains MATLAB code used for the quantitative EEG analyses, biomarker breakpoint analyses, and statistical analyses reported in the accompanying study.

The analysis code is organised into three main components:

1. **EEG identifiability analyses**
2. **Piecewise linear modelling of biomarker trajectories**
3. **Statistical analyses of EEG, clinical, biomarker, and neuropsychological measures**

The code is intended to document the analytical procedures used in the study and to support reproducibility of the reported analyses.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Software requirements

The analyses were developed and run using:

* MATLAB R2023a
* Statistics and Machine Learning Toolbox

Additional MATLAB functions may be required depending on the specific analysis and plotting routines.

The `bh_fdr.m` function included in this repository implements the Benjamini–Hochberg false-discovery-rate correction used in the statistical analyses.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 1. EEG identifiability analysis

## Overview

The EEG identifiability analysis quantifies the stability and distinctiveness of individual spectral profiles across EEG recordings.

For each EEG recording, resting-state EEG epochs were divided into two sets of equal size. Spectral features were calculated separately for the two sets, producing feature vectors for each recording.

The analysis then:

1. Constructs two spectral feature matrices corresponding to the two epoch sets.
2. Z-scores spectral features across recordings.
3. Calculates Pearson correlations between all pairs of feature vectors.
4. Constructs an EEG-by-EEG similarity matrix.
5. Calculates the self-rank of each recording.
6. Quantifies identity strength relative to a between-subject null distribution.
7. Compares within-recording, within-subject, and between-subject similarity using a linear mixed-effects model.
8. Quantifies longitudinal changes in individual spectral similarity.
9. Estimates one longitudinal slope per participant.
10. Compares slopes between diagnostic groups using non-parametric permutation testing.
11. Estimates uncertainty around group mean slopes using bootstrap resampling.
12. Generates the figures used to visualise spectral identifiability and longitudinal trajectories.

## Spectral feature construction

The analysis uses PSD features extracted from five predefined regions:

* central
* left frontal
* right frontal
* left temporal
* right temporal

The resulting feature vector contains:

**5 regions × 129 frequency features = 645 features per recording.**

The code expects the feature matrices to contain features × recordings.

The two matrices are referred to as:

featureA
featureB


where each column corresponds to an EEG recording.

The random split of epochs into the two sets used to generate `featureA` and `featureB` was performed upstream during EEG processing and is therefore not repeated by the identifiability analysis code.

## Group coding

Diagnostic groups are coded as:

0 = control/asymptomatic
1 = converter
2 = symptomatic


The analysis code explicitly assigns these numerical codes to the corresponding group labels.

## Similarity matrix

Pearson correlations are calculated between every feature vector in set A and every feature vector in set B.

The resulting matrix contains:

* diagonal elements: within-recording similarity
* off-diagonal elements from the same participant: within-subject similarity
* off-diagonal elements from different participants: between-subject similarity

## Self-rank

For each EEG recording, the within-recording correlation is ranked against its correlations with all other recordings.

A lower numerical rank therefore indicates a more distinctive/self-identifying spectral profile.

## Identity strength

Identity strength is calculated relative to the between-subject null distribution:

identity strength = (r − mean_between) / SD_between

where `r` is the relevant correlation coefficient.

This null-standardised measure is used for the comparison of within-recording, within-subject, and between-subject similarity.

## Mixed-effects model

The three similarity conditions are compared using a linear mixed-effects model with condition as a fixed effect and crossed random intercepts for the two members of each EEG pair:

Z ~ Condition + (1|ID_i) + (1|ID_j)

This accounts for the dependence introduced by repeatedly occurring participants in the pairwise similarity data.

## Longitudinal analysis

For participants with multiple EEG recordings, the earliest available recording is treated as baseline.

The baseline point is defined using the within-recording similarity of the two spectral feature halves.

For subsequent recordings, similarity is calculated between the baseline recording and each later recording.

Pearson correlations are Fisher-z transformed using:

atanh(r)

For each participant, a simple linear regression of Fisher-z-transformed similarity against time provides one longitudinal slope estimate.

Participants with fewer than two recordings are not included in the participant-level slope analysis.

## Group comparison

Participant-level slopes are compared between all possible pairs of diagnostic groups:

* control vs converter
* control vs symptomatic
* converter vs symptomatic

Group differences are assessed using a non-parametric permutation test with 10,000 permutations.

## Bootstrap confidence intervals

Uncertainty around the group mean slopes is estimated using 5,000 bootstrap resamples.

The reported confidence intervals are percentile-based 95% bootstrap confidence intervals.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 2. Piecewise linear modelling

## Overview

The piecewise modelling analysis investigates whether biomarker trajectories exhibit a change in slope at a particular time point.

A demo.m script which generates synthetic biomarker data is also provided

The analysis is performed in two stages.

### Stage 1 — patient + control model

Longitudinal biomarker data from patients/converters and controls are modelled together.

The purpose of this analysis is to identify biomarkers showing a significant change in trajectory associated with the disease group.

### Stage 2 — patient-only model

For biomarkers identified in Stage 1, a patient-only piecewise model is fitted to estimate the timing of the breakpoint.

### Stage 3 — bootstrap breakpoint comparison

The estimated patient-only breakpoints are compared between biomarkers using paired subject-level bootstrap resampling.

Only subjects with data available for both biomarkers being compared contribute to a paired comparison.

---

## Breakpoint estimation

A single breakpoint, `tau`, is estimated for each biomarker.

Candidate breakpoints are evaluated over a predefined grid spanning the observed range of time up to the reference time point.

For each candidate breakpoint, the piecewise model is fitted using maximum likelihood.

The breakpoint giving the maximum likelihood is selected as the estimated breakpoint:

tau_hat = breakpoint with maximum log likelihood

The likelihood profile is retained to quantify uncertainty around the estimated breakpoint.

## Profile-likelihood confidence interval

A 95% confidence interval for the breakpoint is obtained from the profile likelihood using the likelihood-ratio criterion:

2 × (LLmax − LL(tau)) ≤ χ²(0.95, 1)

equivalently,

2 × (LLmax − LL(tau)) ≤ 3.84

The confidence interval is determined from the range of candidate breakpoints satisfying this criterion.


## Mixed-effects model

The piecewise models account for repeated observations within participants using participant-level random effects.

The model includes:

* age
* sex
* time (relative to clinical disease onset/last follow-up)
* post-breakpoint time
* group and group-by-post-breakpoint-time interaction where applicable
* participant-level random intercept
* participant-level random post-breakpoint slope

The post-breakpoint time variable is defined relative to the estimated breakpoint.

## Biomarker selection

The patient + control model is used to identify biomarkers for which the group-by-post-breakpoint-time interaction is statistically significant after the prespecified multiple-comparison correction.

The resulting significant biomarker list is then supplied to the patient-only breakpoint analysis.

## Bootstrap breakpoint comparison

For selected biomarkers, breakpoints are compared using a paired subject-level bootstrap.

For each biomarker pair:

1. Participants with observations for both biomarkers are identified.
2. Participants are sampled with replacement.
3. The same sampled participant sequence is used for both biomarkers, preserving the paired structure.
4. The patient-only piecewise model is fitted separately to each bootstrap sample.
5. The difference between the estimated breakpoints is calculated.
6. This procedure is repeated 2,000 times.
7. A percentile-based 95% confidence interval is calculated for the breakpoint difference.

The direction of the breakpoint difference is defined by the order of the biomarker pair:

Δtau = tau_A − tau_B

A confidence interval that excludes zero indicates evidence for a difference in estimated breakpoint timing.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 3. Statistical analysis

## Overview

The statistical analysis code contains analyses of:

* baseline quantitative EEG measures
* qualitative clinical EEG findings
* longitudinal quantitative EEG measures
* associations between EEG measures and plasma/neuropsychological variables

## Baseline quantitative EEG

Baseline quantitative EEG measures are compared between:

1. asymptomatic mutation carriers more than two years from expected onset
2. converters within two years of onset
3. symptomatic carriers

One baseline EEG recording per participant is used.

Group differences are assessed using ANCOVA with:

* diagnostic group as the factor
* age as a covariate
* sex as a covariate

Pairwise group comparisons are performed relative to the prespecified reference group where applicable.

## Multiple-comparison correction

Multiple comparisons for quantitative spectral measures are controlled using the Benjamini–Hochberg false discovery rate procedure.

The implementation is contained in:

bh_fdr.m

The function handles missing p-values and returns adjusted q-values.

The definition of each statistical testing family should follow the prespecified analysis plan and the corresponding Methods description.

## Clinical EEG findings

Clinical EEG reports were reviewed by experienced neurophysiologists.

The qualitative analysis classifies clinical EEG findings according to the predefined normal/abnormal or nonspecific categories used in the study.

The association between diagnostic group and qualitative EEG status is evaluated using logistic regression with a binomial distribution and logit link.

Group contrasts are evaluated relative to the asymptomatic (control) reference group.

## Longitudinal quantitative EEG

All available longitudinal EEG recordings from converters are included.

Longitudinal EEG measures are analysed using linear mixed-effects models.

The model includes:

* time (relative to clinical onset, where clinical onset is time 0)
* age
* sex
* participant-level random intercept

The general model is:

EEG metric ~ Time + Age + Sex + (1|Subject)

The random intercept accounts for repeated measurements within participants.

## EEG–plasma/neuropsychological associations

Associations between EEG measures and plasma or neuropsychological variables are analysed using linear mixed-effects models.

The model includes:

* biomarker/neuropsychological measure as predictor
* age
* sex
* participant-level random intercept

All available paired observations are included subject to the prespecified inclusion criteria.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 4. Data availability and confidentiality

Raw EEG recordings and participant-level clinical data are not included in this repository because of confidentiality and data-governance restrictions.

The analysis scripts therefore operate on derived input data that are not publicly distributed.

Where appropriate, synthetic or simulated data may be provided to demonstrate the expected input structure and to allow users to test the analysis functions without access to the restricted study data.

The absence of participant-level data from this repository does not indicate that the analyses require proprietary software or inaccessible analytical procedures; the statistical procedures are implemented using standard MATLAB functionality and the supplied analysis code.

Users attempting to reproduce the numerical results reported in the manuscript will require access to the corresponding study data under the applicable data-sharing and ethical governance procedures.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 5. Citation

If this repository is used or adapted, please cite the accompanying publication.

The EEG spectral identifiability analysis follows previously published correlation-based spectral fingerprinting/identifiability methodology as described in the manuscript.

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

# 6. Contact

For questions regarding the analysis code, input-data structure, or reproducibility of the reported analyses, please refer to the corresponding author or repository maintainers.
