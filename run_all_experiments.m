%========================================
% Source: run_all_experiments.m
% Author: Keith Bush
% Date: Aug 27, 2012
% 
% Revised: February 21, 2014
%
% Revised: April 8, 2014 (rewrote to sectionalize code)
%
% Purpose: Examine the use of resampled deconvolution to strip
% uncertain neural events and reconvolve a "confident BOLD signal".
% This code examines base deconvolution filtering against NEV
% filtered variants in terms of functional connectivity analysis.
%
%========================================
clear all; close all;

%===================================
%Initialize the parameter structure
%===================================
addpath(genpath('/home/kabush/lib/spm8'));
addpath(genpath('./anev_src'));
addpath(genpath('./deconvolve_src'));
addpath(genpath('./utility_src'));

%=====================================
% Initialize anev parameter structurecle
%=====================================
anev_params = struct('savefilename', ([]), 'TS', ([]), 'FG', ([]), 'FO', ([]),...
                     'FA', ([]), 'HRF_d', ([]),'SNphys', ([]), 'SNscan', ([]),...
                     'hrfOn', ([]), 'noiseOn', ([]), 'pruneOn', ([]), ...
                     'percentOn', ([]), 'N', ([]), 'ANmu', ([]), ...
                     'ANvar', ([]), 'ANrot', ([]), 'maxLag', ([]), ...
                     'LAG', ([]));

anev_params.TS = 200;                    % Observed trace length (number of points)
anev_params.FG = 1.0;                    % Frequency generated (Hz)
anev_params.FO = 1.0;                    % Frequency observed (Hz)
anev_params.FA = 0.05;                   % Fraction of generation time-steps containing
                                         % neural event
anev_params.HRF_d = 6;                   % HRF time-to-peak
                                         % parameter (HRF_d=6 is
                                         % default, BLN model equivalent is HRF_d = 4)
anev_params.SNphys = 6.0;                % Physiological signal-to-noise ratio of BOLD
anev_params.SNscan = 9.0;                % Scanner signal-to-noise (defined as power ratio)
anev_params.ARphys = 0.75;               % coef. auto-regression (1-step)
anev_params.ARscan = 0.0;                % coef. auto-regression (1-step).
anev_params.hrfOn =  true;               % Is the BOLD generated from HRF? (o/w default
                                         % random Balloon model is used
anev_params.noiseOn = false;             % Is the noise model applied tot he observation?
anev_params.pruneOn = false;             % Is pruning applied to the observation?
anev_params.percentOn = true;            % Is the signal normalized to
                                         % percent signal change?
anev_params.NTrans = 3;                  % Number of kernel widths considered
                                         % transient signal                                                         
%Default 3-node functional network
anev_params.N = 3;                                         % Number of nodes
anev_params.ANmu = [-4 4 -4 -4 -4 4 -4 -4 -4];             % Mean connection strength
anev_params.ANvar = [0 0 0 0 0 0 0 0 0];                   % Var connection strength
anev_params.ANrot = RandOrthMat(length(anev_params.ANmu)); % Random connection rotation
anev_params.maxLag = 0.050;                                % Maximum lag (s)
anev_params.LAG = anev_params.maxLag*ones(anev_params.N);  % Distribution of lag
                                                           
%==============================================
% Initialize deconvolution parameter structure
%==============================================
dcv_params = struct('hrf_lr', ([]),'epsilon', ([]), 'beta', ([]), 'Nresample', ([]));
dcv_params.lr = 0.01;
dcv_params.epsilon = 0.005;
dcv_params.beta = 40;        % ****Taken from Gokce's Results****
dcv_params.Nresample = 30;

%==============================================
% Initialize classify parameter structure
%==============================================
class_params.prob_step = 0.01;
class_params.prob_min = 0.0;
class_params.prob_max = 0.5;

class_params.std_step = 0.01;
class_params.std_min = 0.0;
class_params.std_max = 1.0;

%==================================
% Initialize experiment parameters
%==================================
Nsamples = 30; 

%================================
% Seed random number explicitely

%================================
%rng shuffle
% stream = RandStream('mt19937ar','Seed',sum(100*clock));
% RandStream.setDefaultStream(stream);

%============================================
% Experiments 1
% Examine full confound w/o misspecification
% 3-node (1->2, 2->3)
%=============================================
anev_params.N = 1;
anev_params.TS = 200;  
anev_params.FG = 20.0;
anev_params.FO = 1.0;  
anev_params.noiseOn = true;
anev_params.pruneOn = true;
anev_params.HRF_d = 6;
anev_params.hrfOn = true;

%Net params (not relevant)
anev_params.ANmu = [-20];
anev_params.ANvar = [0];
anev_params.maxLag = 0.050;                                
anev_params.ANrot = RandOrthMat(length(anev_params.ANmu)); 
anev_params.LAG = anev_params.maxLag*ones(anev_params.N);  

knownHRF = false;   %IE: use HRF_d=6 
HRF_d_mu = 6;  
HRF_d_sigma = 0;    %*********** NO CONFOUND ************

%%Need to control deconvolution performance for ARphys
% ARphys_seq = 0.0:0.05:1.0;
ARphys_seq = 0.75;
for i=1:1 %1:numel(ARphys_seq)

    experiment_id = i;

    %%Set auto-regression characteristic
    anev_params.ARphys = ARphys_seq(i); %0.75

    %%Step 1: Generate the data
    run_generate_rand_hrf(experiment_id, Nsamples, anev_params, ...
                          dcv_params, class_params,HRF_d_mu,HRF_d_sigma);

    %%Step 2: Deconvolve the data into neural events
    run_deconvolve(experiment_id,Nsamples,knownHRF);
    
    %%Step 3: Classify the results
    [results] = run_classify(experiment_id,Nsamples);
    save(['./data/results/experiment',num2str(experiment_id),'.mat'],'results');

end
