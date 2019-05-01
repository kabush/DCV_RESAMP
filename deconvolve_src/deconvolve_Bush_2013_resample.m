% Author:      Keith Bush, PhD
% Institution: University of Arkansas at Little Rock
% Date:        Aug. 12, 2013

function [result] = deconvolve_Bush_2013_resample(BLDobs,kernel,nev_lr,epsilon,beta,Nresample)
%-Description
% This function deconvolves the BOLD signal using Bush 2011 method
%
%-Inputs
% BLDobs - observed BOLD signal
% kernel - assumed kernel of the BOLD signal
% nev_lr - learning rate for the assignment of neural events
% epsilon - relative error change (termination condition)
% beta - slope of the sigmoid transfer function
% Nresample - 
%
%-Outputs
% result - struct of important data
%    NEVest  - the base neural event estimate
%    NEVmean - the mean neural event estimate
%    NEVstd - the std dev. of the neural event estimate
%    NEVcupp - the mean (upper limit 95% ci)
%    NEVclow - the mean (lower limit 95% ci)
%    BLDmean - the mean of the BOLD estimate
%    BLDstd - the std dev. of the BOLD estimate
%    BLDcupp - the mean BOLD (upper limit 95% ci)
%    BLDclow - the mean BOLD (lower limit 95% ci)
    
    Kobs = numel(kernel);
    
    %Scale observe via z-scoring
    BLDobs_scale = zscore(BLDobs);

    %Deconvolve observed BOLD
    [NEVdcv] = deconvolve_Bush_2013(BLDobs,kernel,nev_lr,epsilon,beta);
    
    %Reconvolve estimated true BOLD
    [BLDdcv_full] = convolve_dcv_hrf(NEVdcv',kernel);
    BLDdcv_low = BLDdcv_full(Kobs:end);

    %Convert to a z-score representation
    BLDdcv_scale = zscore(BLDdcv_low);
  
    %==============================
    % Apply the resampling method
    %==============================
    NEVvariants = zeros(Nresample,numel(NEVdcv));
    BLDvariants = zeros(Nresample,numel(BLDobs_scale));

    NEVvariants(1,:) = NEVdcv;
    BLDvariants(1,:) = BLDdcv_scale;

    for z=2:Nresample
        
        %Compute residual
        residuals = BLDobs_scale - BLDdcv_scale;

        %Randomize the residual order
        RND_residuals = residuals(randi(numel(residuals),1,numel(residuals)));
        
        %Reapply the residual to the filtered BLDobs
        RNDobs = BLDdcv_scale + RND_residuals;         

        %Deconvolve the variant
        [NEVresidual] = deconvolve_Bush_2013(RNDobs,kernel,nev_lr,epsilon,beta);
        
        %Store encoding of this variant
        NEVvariants(z,:) = NEVresidual';
        
        %Reconvolve to form the filtered BLD of this variant
        BLDdcv_full = convolve_dcv_hrf(NEVresidual',kernel);
        BLDdcv_low = BLDdcv_full(Kobs:end);
        BLDdcv_res = zscore(BLDdcv_low);
        BLDvariants(z,:) = BLDdcv_res;

    end

    %========================================================
    % Compute Performance for regular versus precision-guided
    %========================================================
    
    %%Compute distribution over NEVvariants
    result.NEVest = NEVdcv;
    result.NEVmean = mean(NEVvariants);
    result.NEVstd = std(NEVvariants);
    
    result.NEVcupp = 0*result.NEVest;
    result.NEVclow = 0*result.NEVest;
    
    for i=1:numel(result.NEVcupp)
       [H P CI STATS] = ttest(NEVvariants(:,i));
       result.NEVclow(i) = CI(1);
       result.NEVcupp(i) = CI(2);
    end
    
    %%Compute distribution over BLDvariants
    result.BLDmean = mean(BLDvariants);
    result.BLDstd = std(BLDvariants);
    result.BLDcupp = 0*result.BLDmean;
    result.BLDclow = 0*result.BLDmean;
    
    for i=1:numel(result.BLDcupp)
       [H P CI STATS] = ttest(BLDvariants(:,i));
       result.BLDclow(i) = CI(1);
       result.BLDcupp(i) = CI(2);
    end

    

end
