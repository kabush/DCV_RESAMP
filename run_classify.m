function [results] = run_classify(experiment_id,Nsamples)

%Allocate storage for analysis (**** THIS IS A HACK!!!!, class and
%anev params should be able to vary trial to trial!!!!!!****)
    class_save = load(['./data/class_data/class_save_',num2str(experiment_id),'_',num2str(1),'.txt']);
    class_params.prob_step = class_save(1);
    class_params.prob_min = class_save(2);
    class_params.prob_max = class_save(3);
    class_params.std_step = class_save(4);
    class_params.std_min = class_save(5);
    class_params.std_max = class_save(6);
    
    %%Reconstruct anev_params
    anev_save = load(['./data/gen_data/anev_save_',num2str(experiment_id),'_',num2str(1),'.txt']);
    anev_params.N = anev_save(2);
    anev_params.TS = anev_save(3);
    anev_params.FG = anev_save(4);
    anev_params.FO = anev_save(5);
    anev_params.HRF_d = anev_save(6);

    %%Std seq values
    std_seq = class_params.std_min:class_params.std_step: ...
              class_params.std_max;
    auc_std_results = zeros(Nsamples,anev_params.N,numel(std_seq));
    frac_std_results = zeros(Nsamples,anev_params.N,numel(std_seq));

    %%Prob seq values
    prob_seq = class_params.prob_min:class_params.prob_step: ...
              class_params.prob_max;
    auc_ci_results = zeros(Nsamples,anev_params.N,numel(prob_seq));    
    frac_ci_results = zeros(Nsamples,anev_params.N,numel(prob_seq));
        
    %%================================
    %% Perform all the classifications
    %%================================
    for sample=1:Nsamples

        ['Classify: ',num2str(sample)]
        
        %%Reconstruct anev_params
        anev_save = load(['./data/gen_data/anev_save_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        anev_params.N = anev_save(2);
        anev_params.TS = anev_save(3);
        anev_params.FG = anev_save(4);
        anev_params.FO = anev_save(5);
        anev_params.HRF_d = anev_save(6);
  
        %%DEBUG
        %% class_save = load(['./data/class_data/class_save_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        %% class_params.prob_step = class_save(1);
        %% class_params.prob_min = class_save(2);
        %% class_params.prob_max = class_save(3);
        %% class_params.std_step = class_save(4);
        %% class_params.std_min = class_save(5);
        %% class_params.std_max = class_save(6);
        
        %from simulation
        inputs.FN = load(['./data/gen_data/FN_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.BLDobs = load(['./data/gen_data/BLDobs_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.BLDobs_adj = load(['./data/gen_data/BLDobs_adj_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.NEVgen_prn = load(['./data/gen_data/NEVgen_prn_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        
        %from deconvolution
        inputs.NEVest = load(['./data/dcv_data/NEVest_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.NEVmean = load(['./data/dcv_data/NEVmean_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.NEVstd = load(['./data/dcv_data/NEVstd_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.NEVcupp = load(['./data/dcv_data/NEVcupp_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        inputs.NEVclow = load(['./data/dcv_data/NEVclow_',num2str(experiment_id),'_',num2str(sample),'.txt']);
        
        %Classify
        [results_MatClassify] = run_MatClassify(inputs,anev_params,class_params);
        
        %%
        %Extract important data
        %%
        auc_std_results(sample,:,:) = results_MatClassify.auc_std_conf;
        frac_std_results(sample,:,:) = results_MatClassify.frac_std_conf;

        auc_ci_results(sample,:,:) = results_MatClassify.auc_ci_conf;
        frac_ci_results(sample,:,:) = results_MatClassify.frac_ci_conf;
        
    end
    
    %%================================
    %% Perform all the analysis
    %%================================    
    auc_std_mean = zeros(anev_params.N,numel(std_seq));
    auc_std_std = zeros(anev_params.N,numel(std_seq));
    auc_std_cupp = zeros(anev_params.N,numel(std_seq));
    auc_std_clow = zeros(anev_params.N,numel(std_seq));
    
    frac_std_mean = zeros(anev_params.N,numel(std_seq));
    frac_std_std = zeros(anev_params.N,numel(std_seq));
    frac_std_cupp = zeros(anev_params.N,numel(std_seq));
    frac_std_clow = zeros(anev_params.N,numel(std_seq));

    auc_ci_mean = zeros(anev_params.N,numel(prob_seq));
    auc_ci_std = zeros(anev_params.N,numel(prob_seq));
    auc_ci_cupp = zeros(anev_params.N,numel(prob_seq));
    auc_ci_clow = zeros(anev_params.N,numel(prob_seq));
    
    frac_ci_mean = zeros(anev_params.N,numel(prob_seq));
    frac_ci_std = zeros(anev_params.N,numel(prob_seq));
    frac_ci_cupp = zeros(anev_params.N,numel(prob_seq));
    frac_ci_clow = zeros(anev_params.N,numel(prob_seq));

    
    for n=1:anev_params.N

        %% ==================================
        %% Summarize AUC for std dev analysis
        %% ==================================
        
        auc_std_mean(n,:) = mean(reshape(auc_std_results(:,n,:), ...
                                         Nsamples,numel(std_seq)));
        
        frac_std_mean(n,:) = mean(reshape(frac_std_results(:,n,:), ...
                                          Nsamples,numel(std_seq)));
        
        for i=1:numel(std_seq)

            %%Compute std dev.
            auc_std_std(n,i) = std(reshape(auc_std_results(:,n,i), ...
                                           Nsamples,1));
            
            %%Compute 95% interval
            [H,P,CI] = ttest(reshape(auc_std_results(:,n,i),Nsamples,1));
            auc_std_clow(n,i) = CI(1);
            auc_std_cupp(n,i) = CI(2);
            
            %%
            %% Summarize Fraction of NEV 
            %%

            %%Compute std dev.
            frac_std_std(n,i) = std(reshape(frac_std_results(:,n,i), ...
                                            Nsamples,1));
            
            %%Compute 95% interval
            [H,P,CI] = ttest(reshape(frac_std_results(:,n,i),Nsamples,1));
            frac_std_clow(n,i) = CI(1);
            frac_std_cupp(n,i) = CI(2);

        end
    
    
        %% ==================================
        %% Summarize AUC for ci analysis
        %% ==================================
        
        auc_ci_mean(n,:) = mean(reshape(auc_ci_results(:,n,:), ...
                                         Nsamples,numel(prob_seq)));

        frac_ci_mean(n,:) = mean(reshape(frac_ci_results(:,n,:), ...
                                         Nsamples,numel(prob_seq)));
    
        for i=1:numel(prob_seq)

            %%Compute std dev.
            auc_ci_std(n,i) = std(reshape(auc_ci_results(:,n,i), ...
                                           Nsamples,1));
            
            %%Compute 95% interval
            [H,P,CI] = ttest(reshape(auc_ci_results(:,n,i),Nsamples,1));
            auc_ci_clow(n,i) = CI(1);
            auc_ci_cupp(n,i) = CI(2);
            
            %%
            %% Summarize Fraction of NEV 
            %%

            %%Compute std dev.
            frac_ci_std(n,i) = std(reshape(frac_ci_results(:,n,i), ...
                                            Nsamples,1));
            
            %%Compute 95% interval
            [H,P,CI] = ttest(reshape(frac_ci_results(:,n,i),Nsamples,1));
            frac_ci_clow(n,i) = CI(1);
            frac_ci_cupp(n,i) = CI(2);

        end
    
    end
    
    %%================================
    %% Gather results
    %%================================            
    results.auc_std_mean = auc_std_mean;
    results.auc_std_std = auc_std_std;
    results.auc_std_cupp = auc_std_cupp;
    results.auc_std_clow = auc_std_clow;
    results.frac_std_mean = frac_std_mean;
    results.frac_std_std = frac_std_std;
    results.frac_std_cupp = frac_std_cupp;
    results.frac_std_clow = frac_std_clow;

    results.auc_ci_mean = auc_ci_mean;
    results.auc_ci_std = auc_ci_std;
    results.auc_ci_cupp = auc_ci_cupp;
    results.auc_ci_clow = auc_ci_clow;
    results.frac_ci_mean = frac_ci_mean;
    results.frac_ci_std = frac_ci_std;
    results.frac_ci_cupp = frac_ci_cupp;
    results.frac_ci_clow = frac_ci_clow;
    
end
