function [results] = run_MatGenerate(anev_params,dcv_params)

%Grab the ANEV params (base signal properties)
TS = anev_params.TS;
FG = anev_params.FG;
FO = anev_params.FO;
FA = anev_params.FA;
HRF_d = anev_params.HRF_d;
SNphys = anev_params.SNphys;
SNscan = anev_params.SNscan;
ARphys = anev_params.ARphys;
ARscan = anev_params.ARscan;

%ANEV params (technical variations)
NTrans = anev_params.NTrans;
hrfOn = anev_params.hrfOn;
noiseOn = anev_params.noiseOn;
pruneOn = anev_params.pruneOn;
percentOn = anev_params.percentOn;
                                         
%Grab the 3-node functional network
N = anev_params.N;
ANmu = anev_params.ANmu;
ANvar = anev_params.ANvar;
ANrot = anev_params.ANrot;
maxLag = anev_params.maxLag;
LAG = anev_params.LAG;

%Build generation kernel
KERgen = spm_advanced_hrf(1/FG,HRF_d);
Kgen = numel(KERgen);
 
%===============================================
%Sample the specific functional network
%===============================================
FS = sample_network(ANmu,ANvar,ANrot);
FN = reshape(FS,N,N)';
B = repmat(FA,N,1);

%===============================================
%Simulate the network to generate neural events
%===============================================
NEVgen = generate_anev_roi_model(FN,LAG,TS,NTrans*Kgen-1,FG,B);

%Allocate tmp storage
NEValt = 0*NEVgen;
BLDgen = 0*NEVgen;
SCKSgen = zeros(7,N);

%Convolve ANEV with the HRF Model
if hrfOn 
    %%Convolve all at once (using canonical HRF)
    BLDgen = convolve_anev_roi_hrf(NEVgen,FG,HRF_d);
else
    %%Convolve each node individually (using Balloon model)
    for i = 1:N
        [NEValt(i,:),BLDgen(i,:),SCKS] = gen_BOLD_from_NEV_Balloon(NEVgen(i,:),FG,FO);
        SCKSgen(:,i) = SCKS.M(1).pE;
    end
end

%Observe the ANEV Model (introduce confounds)
[NEVgen_prn,BLDgen_prn,BLDobs] = observe_roi(NEVgen,BLDgen,FG, ...
                                             FO,SNphys,SNscan, ...
                                             ARphys,ARscan, ...
                                             NTrans*Kgen-1, ...
                                             noiseOn,pruneOn,percentOn);


%Observe the ANEV Model (without confound for analysis)
[NEVgen_prn,BLDgen_prn_adj,BLDobs_adj] = observe_roi(NEVgen, ...
                                                  BLDgen,FG,FO, ...
                                                  SNphys, ...
                                                  SNscan, ...
                                                  ARphys, ...
                                                  ARscan, ...
                                                  NTrans*Kgen- ...
                                                  1,false,pruneOn,percentOn);

%% RESULTS: NEVgen, BLDgen, NEVgen_prn, BLDobs, BLDobs_adj, HRF_d 
results.FN = FN;
results.NEVgen = NEVgen;
results.BLDgen = BLDgen;
results.NEVgen_prn = NEVgen_prn;
results.BLDgen_prn = BLDgen_prn;
results.BLDobs = BLDobs;
results.BLDobs_adj = BLDobs_adj;
results.HRF_d = HRF_d;