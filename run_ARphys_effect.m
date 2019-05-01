%========================================
% Source: run_ARphys_effect.m
% Author: Keith Bush
% Date: 
%========================================
clear all; close all;
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
Nexp = 21
Nbins = 51

auc_results = zeros(Nexp,Nbins);
frac_results = zeros(Nexp,Nbins);

for i=1:Nexp

    experiment_id = i;
    
    load(['./data/results/experiment',num2str(experiment_id),'.mat'],'results');
    
    %%Analysis
    auc_results(i,:) = results.auc_ci_mean;
    frac_results(i,:) = results.frac_ci_mean;

end

Nadj = Nbins-2
slopes = zeros(Nbins-2,1)
ys = zeros(Nbins-2,1)

for i=1:(Nadj)

    [b_fit stats_fit] = robustfit([0.0:0.05:1.0],auc_results(:,i))
    slopes(i) = b_fit(2)
    ys(i) = b_fit(1)
 
end

%Contruct fit
indep = (1:Nadj)/(2*Nadj)
[b_meta_fit stats_meta_fit] = robustfit(indep,slopes(1:Nadj))
fit_x = (0.0:0.01:0.48)
fit_line = b_meta_fit(1)+b_meta_fit(2)*fit_x%indep

%Plot fit
h=figure(1)
plot(fit_x,slopes(1:Nadj))
hold on;
plot(fit_x,fit_line,'r')
xlabel('\delta')
ylabel('Slope of AUC vs ARphys')
r2=rsquare(slopes,fit_line')
text(0.3,-0.004,['r^2 = ',num2str(r2)])

%Write out plot to file
print(h,'-deps','figAppB.eps')
