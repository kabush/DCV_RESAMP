function [results] = run_MatClassify(inputs,anev_params,class_params)

    %Simulation Params
    N = anev_params.N;
    TS = anev_params.TS;
    FG = anev_params.FG;
    FO = anev_params.FO;
    HRF_d = anev_params.HRF_d;

    FN = inputs.FN;
    
    %Build observation kernel
    KERobs = spm_advanced_hrf(1/FO,HRF_d);
    Kobs = numel(KERobs);
    
    %Neural analysis storage
    NEVest_low = inputs.NEVest;
    NEVmean_low = inputs.NEVmean;
    NEVstd_low = inputs.NEVstd;
    NEVcupp_low = inputs.NEVcupp;
    NEVclow_low = inputs.NEVclow;

    %NEV
    NEVgen_prn = inputs.NEVgen_prn;
    
    %BOLD
    BLDobs = inputs.BLDobs;
    BLDobs_adj = inputs.BLDobs_adj;

    %%Allocate storage for full fidelity signals
    tmp = expand_encoding(NEVest_low(1,Kobs:end),FG,FO);
    NEVest_rec = zeros(N,numel(tmp));
    NEVmean_rec = zeros(N,numel(tmp));
    NEVstd_rec = zeros(N,numel(tmp));
    NEVclow_rec = zeros(N,numel(tmp));
    NEVcupp_rec = zeros(N,numel(tmp));

    %AUC Classification
    prob_seq = class_params.prob_min:class_params.prob_step:class_params.prob_max;
    Nseq = numel(prob_seq); 

    auc_ci_conf = zeros(N,Nseq);
    frac_ci_conf = zeros(N,Nseq);

    std_seq = class_params.std_min:class_params.std_step:class_params.std_max;
    Nstd_seq = numel(std_seq);
    auc_std_conf = zeros(N,Nstd_seq);
    frac_std_conf = zeros(N,Nstd_seq);
 
    %Storage for Error Calculations (PASS OUT)
    errorBLD_wrt_TBLD = zeros(Nseq);
    errorBLD_wrt_FN = zeros(Nseq);
    errorDCV_wrt_TBLD = zeros(Nseq);
    errorDCV_wrt_FN = zeros(Nseq);
    frac_data = zeros(N,Nseq);

    %===============================================
    % Reconstruct full fidelity signals
    %===============================================
    for i = 1:N
        NEVest_rec(i,:) = expand_encoding(NEVest_low(i,Kobs:end),FG,FO);
        NEVmean_rec(i,:) = expand_encoding(NEVmean_low(i,Kobs:end),FG,FO);
        NEVstd_rec(i,:) = expand_encoding(NEVstd_low(i,Kobs:end), ...
                                          FG,FO);
        NEVclow_rec(i,:) = expand_encoding(NEVclow_low(i,Kobs:end), ...
                                           FG,FO);
        NEVcupp_rec(i,:) = expand_encoding(NEVcupp_low(i,Kobs:end),FG,FO);
    end
    
    %=============================================
    % Compute classification accuracy (base/mean)
    %=============================================
    auc_base = zeros(N,1);
    auc_mean = zeros(N,1);
    
    for i = 1:N
        
        %%compute auc of base
        [ROCx,ROCy] = compute_roc(NEVest_rec(i,:),NEVgen_prn(i,:));
        auc_base(i) = trapz(fliplr(ROCx),fliplr(ROCy));
        
        %%compute auc of resampled mean
        [ROCx,ROCy] = compute_roc(NEVmean_rec(i,:),NEVgen_prn(i,:));
        auc_mean (i) = trapz(fliplr(ROCx),fliplr(ROCy));
        
    end
    
    auc_base
    auc_mean
    
    %===========================================================
    % Assess classification performance (resampling) by variance
    %===========================================================    
    
    for i = 1:N
        
        if i==1
            figure(101)
            plot(NEVmean_rec(i,:),'LineWidth',2);
            hold on;
            plot(NEVmean_rec(i,:)+3*NEVstd_rec(i,:));
            plot(NEVmean_rec(i,:)-3*NEVstd_rec(i,:));
            plot(NEVest_rec(i,:),'r');
            plot(NEVgen_prn(i,:)*1.0,'c');
            hold off;
            drawnow;
        end
        
        NEVstd_conf = NEVstd_rec(i,:)/max(NEVstd_rec(i,:));
        
        for j = 1:Nstd_seq
            
            %%Select neural events by confidence interval
            keep_ids = find(NEVstd_conf<std_seq(j));
                            
            %%Compute fraction of elements at threshold
            frac_std_conf(i,j) = numel(keep_ids)/numel(NEVstd_conf);
            
            %%Compute auc of elements at threshold
            [ROCx,ROCy] = compute_roc(NEVest_rec(i,keep_ids),NEVgen_prn(i,keep_ids));
            auc_std_conf(i,j) = trapz(fliplr(ROCx),fliplr(ROCy));
            
        end
    end

    %auc_std_conf
    
    %=======================================================
    % Assess classification performance (resampling) by conf
    %=======================================================
    
    for i = 1:N
        
        %% figure(100)
        %% plot(NEVclow_rec(i,:));
        %% hold on;
        %% plot(NEVcupp_rec(i,:));
        %% plot(NEVest_rec(i,:),'r');
        %% plot(NEVgen_prn(i,:)*1.0,'c');
        %%
        %% hold off;
        %% pause;
        
        for j = 1:Nseq
        
            %%Select neural events above upper conf. interval
            keep_ids = find(NEVclow_rec(i,:)>=0.5+prob_seq(j));
            
            %%Select neural events below lower conf. interval
            keep_ids = [keep_ids,find(NEVcupp_rec(i,:)<=0.5- ...
                                      prob_seq(j))];
            
            %%Compute fraction of elements at threshold
            frac_ci_conf(i,j) = numel(keep_ids)/numel(NEVest_rec(i,:));
            
            %%Compute auc of elements at threshold
            [ROCx,ROCy] = compute_roc(NEVest_rec(i,keep_ids),NEVgen_prn(i,keep_ids));
            auc_ci_conf(i,j) = trapz(fliplr(ROCx),fliplr(ROCy));

            %% %%Select neural events above upper conf. interval
            %% keep_ids = find((NEVmean_rec(i,:)-3*NEVstd_rec(i,:))>=0.5+prob_seq(j));
            %% 
            %% %%Select neural events below lower conf. interval
            %% keep_ids = [keep_ids,find((NEVmean_rec(i,:)+3*NEVstd_rec(i,:))<=0.5- ...
            %%                           prob_seq(j))];
            %%
            %% %%Compute fraction of elements at threshold
            %% frac_ci_conf(i,j) = numel(keep_ids)/numel(NEVest_rec(i,:));
            %% 
            %% %%Compute auc of elements at threshold
            %% [ROCx,ROCy] = compute_roc(NEVest_rec(i,keep_ids),NEVgen_prn(i,keep_ids));
            %% auc_ci_conf(i,j) = trapz(fliplr(ROCx),fliplr(ROCy));

            
        end
    end
    
    %auc_ci_conf

    %===============================================
    % Filter Neural Events by confidence
    %===============================================
    
    if N>1
        
        %%All filtered indices
        NEVids_low = zeros(N,Nseq,TS*FO+(Kobs-1));
        
        for i=1:N
            
            for j=1:Nseq
                
                %%Select neural events above upper conf. interval
                keep_ids = find(NEVclow_low(i,:)>=0.5+prob_seq(j));
                
                %%Select neural events below lower conf. interval
                keep_ids = [keep_ids,find(NEVcupp_low(i,:)<=0.5- ...
                                          prob_seq(j))];
                
                NEVids_low(i,j,keep_ids) = 1;
                
            end
        end
        
        %===============================================
        % Assess func. connect. id performance
        %===============================================
        True_corr = corr(BLDobs_adj');
        Obs_corr = corr(BLDobs');
        
        Nelements = numel(NEVids_low(1,1,:));
        
        for k=1:Nseq
            
            count = 1;
            this_vec = zeros(1,(N^2-N)/2);
            trueBLD_vec = zeros(1,(N^2-N)/2);
            trueFN_vec = zeros(1,(N^2-N)/2);
            obsBLD_vec = zeros(1,(N^2-N)/2);
            
            for i=1:N
                
                for j=1:N
                    
                    if j>i
                        
                        %%Find ids that intersect at the confidence level
                        conf_ids = intersect(find(NEVids_low(i,k,:)==1),find(NEVids_low(j,k,:)==1));
                        
                        %%Make sure only part of observed time-series
                        conf_ids = intersect(conf_ids,Kobs:Nelements);
                        
                        if(numel(conf_ids)>0)
                            this_vec(count) = corr(NEVest_low(i,conf_ids)',NEVest_low(j,conf_ids)');
                        end
                        
                        %%Grab measures of comparison
                        trueBLD_vec(count) = True_corr(i,j);
                        trueFN_vec(count) = FN(i,j);
                        obsBLD_vec(count) = Obs_corr(i,j);
                        
                        %%Grab fraction of data
                        frac_data(count,k) = numel(conf_ids)/Nelements;
                        
                        %%Update vector id
                        count = count + 1;
                        
                    end
                    
                end        
            end
            
            %Compare with True BOLD func. connectivity
            errorBLD_wrt_TBLD(k) = sqrt(sum((obsBLD_vec-trueBLD_vec).^2));
            errorDCV_wrt_TBLD(k) = sqrt(sum((this_vec-trueBLD_vec).^2));
            
            %Compare with generating model
            errorBLD_wrt_FN(k) = sqrt(sum((obsBLD_vec-trueFN_vec).^2));
            errorDCV_wrt_FN(k) = sqrt(sum((this_vec-trueFN_vec).^2));
            
        end
        
    end

    %Output Results
    results.errorBLD_wrt_TBLD = errorBLD_wrt_TBLD;
    results.errorBLD_wrt_FN = errorBLD_wrt_FN;
    results.errorDCV_wrt_TBLD = errorDCV_wrt_TBLD;
    results.errorDCV_wrt_FN = errorDCV_wrt_FN;
    
    results.frac_data = frac_data;
    results.frac_std_conf = frac_std_conf;
    results.frac_ci_conf = frac_ci_conf;

    results.auc_base = auc_base;
    results.auc_mean = auc_mean;
    results.auc_std_conf = auc_std_conf;
    results.auc_ci_conf = auc_ci_conf;

end
