%========================================
% Source: run_all_analysis.m
% Author: Keith Bush
% Date: Aug 4, 2014
% 
% Purpose: Analysis resampling experiment data.Examine the use of
% resampled deconvolution to strip
%========================================
clear all; close all;

experiment_id = 17;

%===================================
%Initialize the parameter structure
%===================================
addpath(genpath('C:\cygwin\home\kabush\MATLAB\spm8\'));
addpath(genpath('./anev_src'));
addpath(genpath('./deconvolve_src'));
addpath(genpath('./utility_src'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['./data/results/experiment',num2str(experiment_id),'.mat'],'results');

%%Analysis
auc_std_mean = results.auc_std_mean;
auc_std_std = results.auc_std_std;
auc_std_cupp = results.auc_std_cupp;
auc_std_clow = results.auc_std_clow;
frac_std_mean = results.frac_std_mean;
frac_std_std = results.frac_std_std;
frac_std_cupp = results.frac_std_cupp;
frac_std_clow = results.frac_std_clow;

auc_ci_mean = results.auc_ci_mean;
auc_ci_std = results.auc_ci_std;
auc_ci_cupp = results.auc_ci_cupp;
auc_ci_clow = results.auc_ci_clow;
frac_ci_mean = results.frac_ci_mean;
frac_ci_std = results.frac_ci_std;
frac_ci_cupp = results.frac_ci_cupp;
frac_ci_clow = results.frac_ci_clow;

Nstart = 5; %first non-NAN index of auc_std_mean
Nend = 101; %%size of std_seq

%% h = figure(1)
%% plot(Nstart:Nend,auc_std_mean(Nstart:Nend),'LineWidth',2);
%% hold on;
%% plot(Nstart:Nend,auc_std_mean(Nstart:Nend)+auc_std_std(Nstart:Nend));
%% plot(Nstart:Nend,auc_std_mean(Nstart:Nend)-auc_std_std(Nstart:Nend));
%% hold off;
%% print(h,'-dpng','fig1.png');
%% 
%% h=figure(2)
%% plot(Nstart:Nend,frac_std_mean(Nstart:Nend),'Linewidth',2);
%% hold on;
%% plot(Nstart:Nend,frac_std_mean(Nstart:Nend)+frac_std_std(Nstart:Nend));
%% plot(Nstart:Nend,frac_std_mean(Nstart:Nend)-frac_std_std(Nstart:Nend));
%% hold off;
%% print(h,'-dpng','fig2.png');
%% 
%% h = figure(3)
%% plot(Nstart:Nend,auc_std_mean(:,Nstart:Nend),'LineWidth',2);
%% hold on;
%% plot(Nstart:Nend,auc_std_cupp(:,Nstart:Nend));
%% plot(Nstart:Nend,auc_std_clow(:,Nstart:Nend));
%% hold off;
%% print(h,'-dpng','fig3.png');

Nstart =1;
Nend = 50;

h = figure(4)
plot((Nstart:Nend)/100,auc_ci_mean(:,Nstart:Nend),'LineWidth',2);
hold on;
plot((Nstart:Nend)/100,auc_ci_cupp(:,Nstart:Nend));
plot((Nstart:Nend)/100,auc_ci_clow(:,Nstart:Nend));
xlabel('buffer')
ylabel('AUC')
hold off;
print(h,'-dpng','fig4.png');
print(h,'-deps','fig4.eps');

%% h = figure(5)
%% plot(Nstart:Nend,auc_ci_mean(:,Nstart:Nend),'LineWidth',2);
%% hold on;
%% plot(Nstart:Nend,auc_ci_mean(:,Nstart:Nend)+auc_ci_std(:,Nstart:Nend));
%% plot(Nstart:Nend,auc_ci_mean(:,Nstart:Nend)-auc_ci_std(:,Nstart:Nend));
%% xlabel('Threshold')
%% ylabel('AUC')
%% 
%% hold off;
%% print(h,'-dpng','fig5.png');

Nend = 50;

h=figure(6)
plot((Nstart:Nend)/100,frac_ci_mean(Nstart:Nend),'Linewidth',2);
hold on;
plot((Nstart:Nend)/100,frac_ci_cupp(Nstart:Nend));
plot((Nstart:Nend)/100,frac_ci_clow(Nstart:Nend));
xlabel('buffer')
ylabel('Fraction of Known Neural Events')
hold off;
print(h,'-dpng','fig6.png');
print(h,'-deps','fig6.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id = 1
K = 33

NEVest = load(['./data/dcv_data/NEVest_1_',num2str(id),'.txt']);
NEVmean = load(['./data/dcv_data/NEVmean_1_',num2str(id),'.txt']);
NEVclow = load(['./data/dcv_data/NEVclow_1_',num2str(id),'.txt']);
NEVcupp = load(['./data/dcv_data/NEVcupp_1_',num2str(id),'.txt']);
NEVstd = load(['./data/dcv_data/NEVstd_1_',num2str(id),'.txt']);
NEVgen = load(['./data/gen_data/NEVgen_prn_1_',num2str(id),'.txt']);

BLDmean = load(['./data/dcv_data/NEVmean_1_',num2str(id),'.txt']);
BLDclow = load(['./data/dcv_data/NEVclow_1_',num2str(id),'.txt']);
BLDcupp = load(['./data/dcv_data/NEVcupp_1_',num2str(id),'.txt']);

%Prune and expand data to match NEVgen
NEVest = NEVest(K:end);
NEVest_exp = expand_encoding(NEVest,20,1)
NEVmean = NEVmean(K:end);
NEVmean_exp = expand_encoding(NEVmean,20,1)
NEVstd = NEVstd(K:end);

%Variance
NEVstd_exp = expand_encoding(NEVstd,20,1)
NEVsupp_exp = NEVmean_exp+2*NEVstd_exp
NEVslow_exp = NEVmean_exp-2*NEVstd_exp
NEVsupp_exp(find(NEVsupp_exp > 1))=1;
NEVslow_exp(find(NEVslow_exp < 0))=0;

%SEM 
NEVcupp = NEVcupp(K:end);
NEVcupp_exp = expand_encoding(NEVcupp,20,1)
NEVclow = NEVclow(K:end);
NEVclow_exp = expand_encoding(NEVclow,20,1)
NEVcupp_exp(find(NEVcupp_exp > 1))=1;
NEVclow_exp(find(NEVclow_exp < 0))=0;


%% h = figure(7)
%% plot(NEVgen,'c')
%% hold on;
%% plot(NEVest_exp,'color','r')
%% plot(NEVmean_exp,'linewidth',2)
%% plot(NEVsupp_exp)
%% plot(NEVslow_exp)

h = figure(8)
plot(0.995*NEVgen,'c')
hold on;
plot(NEVest_exp,'color','r')
plot(NEVmean_exp,'linewidth',1,'color','b')
plot(NEVcupp_exp,'linewidth',1,'color','b')
plot(NEVclow_exp,'linewidth',1,'color','b')

xlabel('Simulation Steps')
ylabel('Neural Event Probability')
print(h,'-dpng','fig8.png');
print(h,'-deps','fig8.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BLDobs = load(['./data/gen_data/BLDobs_1_',num2str(id),'.txt']);
BLDmean= load(['./data/rsmp_data/BLDmean_1_',num2str(id),'.txt']);
BLDstd = load(['./data/rsmp_data/BLDstd_1_',num2str(id),'.txt']);
BLDcupp = load(['./data/rsmp_data/BLDcupp_1_',num2str(id),'.txt']);
BLDclow = load(['./data/rsmp_data/BLDclow_1_',num2str(id),'.txt']);
BLDsupp = BLDmean+2*BLDstd
BLDslow = BLDmean-2*BLDstd

h = figure(9)
plot(zscore(BLDobs),'r')
hold on;
plot(BLDmean,'linewidth',1,'color','b')
%% plot(BLDcupp,':b')
%% plot(BLDclow,':b')
plot(BLDsupp,'b')
plot(BLDslow,'b')
xlabel('Samples')
ylabel('Observed BOLD')
print(h,'-dpng','fig9.png');
print(h,'-deps','fig9.eps');

%% h = figure(10)
%% plot(zscore(BLDobs),'r')
%% hold on;
%% plot(BLDmean,'linewidth',2,'color','b')
%% plot(BLDsupp,'b')
%% plot(BLDslow,'b')
