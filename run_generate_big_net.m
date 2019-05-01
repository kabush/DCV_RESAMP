%========================================
% Source: run_generate_big_net.m
% Author: Keith Bush
% Date: April 4, 2014
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
addpath(genpath('C:\cygwin\home\kabush\MATLAB\spm8\'));
addpath(genpath('./anev_src'));
addpath(genpath('./deconvolve_src'));
addpath(genpath('./utility_src'));

%%Neural Event Generation
N = 10;
TS = 200;                    % Observed trace length (number of points)
Ttrans = 2*TS;
FG = 20.0;                   % Frequency generated (Hz)
minLag = .010;
maxLag = .050;
Amin = -20;                  %Minimum connection strength
Amax = 4;                    %Maximum connection strength
Bmin = 0.01;                 %Minimum fraction external activity (neural events)
Bmax = 0.01;                 %Maximum fraction external activity (neural events)

%%Observation Generation
FO = 0.5;                    % Frequency observed (Hz)
HRFmu = 6;                   % HRF time-to-peak
                                         % parameter (HRF_d=6 is
                                         % default, BLN model
                                         % equivalent is HRF_d = 4)
HRFsigma = 1;  %%variance over HRFmu
SNphys = 6.0;                % Physiological signal-to-noise ratio of BOLD
SNscan = 9.0;                % Scanner signal-to-noise (defined as power ratio)
ARphys = 0.75;               % coef. auto-regression (1-step)
ARscan = 0.0;                % coef. auto-regression (1-step).

noiseOn = true;             % Is the noise model applied tot he observation?
pruneOn = true;             % Is pruning applied to the observation?
percentOn = true;           % Is the signal normalized to
                            % percent signal change?

%Base System (perfect data, etc)
net = genSimulation(Ttrans+TS,0,N,N,FG,minLag,maxLag,Amin,Amax,Bmin,Bmax);
NEVgen = net.NEVgen;

%% %%Convolve
[BLDgen HRFgen] = convolve_anev_roi_hrf_dist(net.NEVgen,FG,HRFmu,HRFsigma);

%Observe the ANEV Model (truly)
[NEVgen_prn,BLDgen_prn,BLDobs] = observe_roi(NEVgen,BLDgen,FG, ...
                                             FO,SNphys,SNscan, ...
                                             ARphys,ARscan, ...
                                             Ttrans*FG, ...
                                             noiseOn,pruneOn,percentOn);

%Observe the ANEV Model (without confound for analysis)
[NEVgen_prn,BLDgen_prn_adj,BLDobs_adj] = observe_roi(NEVgen, ...
                                                  BLDgen,FG,FO, ...
                                                  SNphys, ...
                                                  SNscan, ...
                                                  ARphys, ...
                                                  ARscan, ...
                                                  Ttrans*FG,false,pruneOn,percentOn);


save(['./data_sim/BLDobs',num2str(N),'.txt'],'-ascii','BLDobs');
save(['./data_sim/BLDobs_adj',num2str(N),'.txt'],'-ascii','BLDobs_adj');
save(['./data_sim/NEVgen_prn',num2str(N),'.txt'],'-ascii','NEVgen_prn');            

figure(1)
subplot(1,2,1)
plot(BLDobs(1,:));
hold on;
for i=2:10
   plot(BLDobs(i,:));
    
end
hold off;

subplot(1,2,2)
plot(BLDobs_adj(1,:));
hold on;
for i=2:10
   plot(BLDobs_adj(i,:));
    
end
hold off;

