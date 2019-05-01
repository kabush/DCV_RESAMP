function run_generate_rand_hrf(experiment_id, Nsamples, anev_params, dcv_params,class_params,HRF_mu,HRF_sigma)
 
    for s=1:Nsamples
        ['Generate: ',num2str(s)]

        [results_MatGenerate,HRF_d_gen] = run_MatGenerateRandHRF(anev_params,dcv_params,HRF_mu,HRF_sigma);
        
        %%Format for output to file
        anev_save = [Nsamples, anev_params.N, anev_params.TS, anev_params.FG, ...
                     anev_params.FO, anev_params.HRF_d, HRF_d_gen];
        dcv_save = [dcv_params.Nresample,dcv_params.lr,dcv_params.epsilon, ...
                    dcv_params.beta];
        
        class_save = [class_params.prob_step, class_params.prob_min, ...
                      class_params.prob_max, class_params.std_step, ...
                      class_params.std_min, class_params.std_max];

        
        FN = results_MatGenerate.FN;
        BLDobs = results_MatGenerate.BLDobs;
        BLDobs_adj = results_MatGenerate.BLDobs_adj;
        NEVgen_prn = results_MatGenerate.NEVgen_prn;
        
        %%Output to file
        save(['./data/gen_data/anev_save_',num2str(experiment_id),'_',num2str(s),'.txt'],'anev_save','-ascii');
        save(['./data/gen_data/dcv_save_',num2str(experiment_id), '_',num2str(s),'.txt'],'dcv_save','-ascii');
        save(['./data/class_data/class_save_',num2str(experiment_id),'_',num2str(s),'.txt'],'class_save','-ascii');
        save(['./data/gen_data/FN_',num2str(experiment_id),'_',num2str(s),'.txt'],'FN','-ascii');
        save(['./data/gen_data/BLDobs_',num2str(experiment_id),'_',num2str(s),'.txt'],'BLDobs','-ascii');
        save(['./data/gen_data/BLDobs_adj_',num2str(experiment_id),'_',num2str(s),'.txt'],'BLDobs_adj','-ascii');
        save(['./data/gen_data/NEVgen_prn_',num2str(experiment_id),'_',num2str(s),'.txt'],'NEVgen_prn','-ascii');
        
        %%Delete from memory
        clear anev_save;
        clear dcv_save;
        clear FN;
        clear BLDobs;
        clear BLDobs_adj;
        clear NEVgen_prn;
        
    end
    
end
