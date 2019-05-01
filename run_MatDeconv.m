function [results] = run_MatDeconv(BLDobs,anev_params,dcv_params,knownHRF)

    %Anev Params
    N = anev_params.N;
    FO = anev_params.FO;
    
    if(knownHRF)
        HRF_d = anev_params.HRF_d_gen;  %%Use exact
    else
        HRF_d = anev_params.HRF_d; %%Use blanket HRF_d (i.e., 6)
    end
    
    %Deconvolution Params
    lr = dcv_params.lr;
    epsilon = dcv_params.epsilon;
    beta = dcv_params.beta;
    Nresample = dcv_params.Nresample;
    
    %Build observation kernel
    KERobs = spm_advanced_hrf(1/FO,HRF_d);

    for i=1:N
        
        %Deconvolve the BOLD (via Bush2013 algorithm + resampling add-on)
        [result] = deconvolve_Bush_2013_resample(BLDobs(i,:), ...
                                                 KERobs,lr,epsilon, ...
                                                 beta, Nresample);

        %Store resampling-deconvolution data for later use
        NEVest_low(i,:) = result.NEVest;
        NEVmean_low(i,:) = result.NEVmean;
        NEVstd_low(i,:) = result.NEVstd;
        NEVcupp_low(i,:) = result.NEVcupp;
        NEVclow_low(i,:) = result.NEVclow;
        
        BLDmean_low(i,:) = result.BLDmean;
        BLDstd_low(i,:) = result.BLDstd;
        BLDcupp_low(i,:) = result.BLDcupp;
        BLDclow_low(i,:) = result.BLDclow;
        
    end
       
    %Neural analysis storage (PASS OUT)
    results.NEVest_low = NEVest_low;
    results.NEVmean_low = NEVmean_low;
    results.NEVstd_low = NEVstd_low;
    results.NEVcupp_low = NEVcupp_low;
    results.NEVclow_low = NEVclow_low;

    %BOLD analysis storage (PASS OUT)
    results.BLDmean_low = BLDmean_low;
    results.BLDstd_low = BLDstd_low;
    results.BLDcupp_low = BLDcupp_low;
    results.BLDclow_low = BLDclow_low;
        
end