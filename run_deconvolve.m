function run_deconvolve(experiment_id,Nsamples,knownHRF)

    ['*********** Deconvolving Data **********']
    
    for s=1:Nsamples

        tStart = tic;
        
        %%Load from file
        BLDobs = load(['./data/gen_data/BLDobs_',num2str(experiment_id),'_',num2str(s),'.txt']);
        anev_save = load(['./data/gen_data/anev_save_',num2str(experiment_id),'_',num2str(s),'.txt']);
        dcv_save = load(['./data/gen_data/dcv_save_',num2str(experiment_id),'_',num2str(s),'.txt']);
        
        %%Reconstruct anev_params
        anev_params.N = anev_save(2);
        anev_params.FO = anev_save(5);
        anev_params.HRF_d = anev_save(6);
        anev_params.HRF_d_gen = anev_save(7);
        
        %Reconstruct dcv_params
        dcv_params.Nresample = dcv_save(1);
        dcv_params.lr = dcv_save(2);
        dcv_params.epsilon = dcv_save(3);
        dcv_params.beta = dcv_save(4);

        %%Run deconvolution
        [results_MatDeconv] = run_MatDeconv(BLDobs,anev_params, ...
                                            dcv_params,knownHRF);
        
        %%Format for output to file
        NEVest = results_MatDeconv.NEVest_low;
        NEVmean = results_MatDeconv.NEVmean_low;
        NEVstd = results_MatDeconv.NEVstd_low;
        NEVcupp = results_MatDeconv.NEVcupp_low;
        NEVclow = results_MatDeconv.NEVclow_low;
        
        %BOLD analysis storage (PASS OUT)
        BLDmean = results_MatDeconv.BLDmean_low;
        BLDstd =  results_MatDeconv.BLDstd_low;
        BLDcupp = results_MatDeconv.BLDcupp_low;
        BLDclow = results_MatDeconv.BLDclow_low;
        
        %%Output NEV to file
        save(['./data/dcv_data/NEVest_',num2str(experiment_id),'_',num2str(s),'.txt'],'NEVest','-ascii');
        save(['./data/dcv_data/NEVmean_',num2str(experiment_id),'_',num2str(s),'.txt'],'NEVmean','-ascii');
        save(['./data/dcv_data/NEVstd_',num2str(experiment_id),'_',num2str(s),'.txt'],'NEVstd','-ascii');
        save(['./data/dcv_data/NEVcupp_',num2str(experiment_id),'_',num2str(s),'.txt'],'NEVcupp','-ascii');
        save(['./data/dcv_data/NEVclow_',num2str(experiment_id),'_',num2str(s),'.txt'],'NEVclow','-ascii');

        %%Output resampled BLD to file
        save(['./data/rsmp_data/BLDmean_',num2str(experiment_id),'_',num2str(s),'.txt'],'BLDmean','-ascii');
        save(['./data/rsmp_data/BLDstd_',num2str(experiment_id),'_',num2str(s),'.txt'],'BLDstd','-ascii');
        save(['./data/rsmp_data/BLDcupp_',num2str(experiment_id),'_',num2str(s),'.txt'],'BLDcupp','-ascii');
        save(['./data/rsmp_data/BLDclow_',num2str(experiment_id),'_',num2str(s),'.txt'],'BLDclow','-ascii');
        
        %Clean-up memory
        clear BLDobs;
        clear NEVest;
        clear NEVmean;
        clear NEVstd;
        clear NEVcupp;
        clear NEVclow;
        clear BLDmean;
        clear BLDstd;
        clear BLDcupp;
        clear BLDclow;
        
        ['Deconvolve: ',num2str(s),' Elapsed: ',num2str(toc(tStart))]
        
    end

end