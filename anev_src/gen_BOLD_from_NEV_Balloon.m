function [NEValt,BLDgen,SCKS] = gen_BOLD_from_NEV_Balloon(NEV,Fgen,Fobs)

% it needs paths to SCKS code
%%    addpath(genpath('./SCKS/'));
    
    %INPUTS
    % NEV - a binaray array of neural events 
    % Fgen - frequency of generation NEV
    % Fobs - frequency that the BOLD will be observed
    % normalize - boolean variable to normalize the BOLD or not
    %OUPUTS
    % SCKS - simulation object containing
    %        M - HDM (see spm8)
    %        Y - array of observations
    %        pU - ????
    %        pP - ????
    %        pH - ????
    
    %set base parameters
    %==========================================================================
    M(1).E.linear = 0;                          % linear model
    M(1).E.s      = 1;                          % smoothness
    M(1).E.dt     = 1/Fgen;                     % integration step for simulation 

    Gkernel = 8; %% Havlicek 2011 parameter (never change)
    
    % level 1
    %------------------------------------------------------------------
    [pE pC] = spm_hdm_priors(1,3);              %
   % pE(6) = 0.02;                               %???? SOURCE ???
    pE(7) = 0.54;                               %???? SOURCE ???
    
    % ========================================
    % ========================================
    % *** randomize first 3 parameters ***
    % ========================================
    % ========================================
    pE([1 2 3]) = pE([1 2 3])+ 1/12*randn(3,1); %%Always randomize
                                                %%Havlicek 2011

    M(1).n  = 4;                                %dim. of state
    M(1).f  = 'spm_fx_hdm';
    M(1).g  = 'spm_gx_hdm';
    M(1).pE = pE;                               % prior expectation
    M(1).V  = exp(6);                           % error precision
    M(1).W  = exp(10);                          % error precision
    
    % level 2
    %------------------------------------------------------------------
    M(2).l  = exp(0);                           % inputs
    M(2).V  = exp(0);                           % with shrinkage priors
    
    M       = spm_DEM_M_set(M);
    
    % free parameters
    %--------------------------------------------------------------------------
    P       = M(1).pE;                                % true parameters
    ip      = [1:length(P)];                          % free parameters
    pE      = spm_vec(P);
    np      = length(pE);
    pE      = spm_unvec(pE,P);
    M(1).pE = pE;
    M(1).pC = pC;
    M(1).ip = ip;
    
    % generate data from SCKS CODE
    %==========================================================================
    N = length(NEV); % length of data
    
    %generate smooth neural inputs
    Ugaussian = (exp(-([1:(2*Gkernel-1)] - Gkernel).^2/(1.^2))); % this is the Gaussian kernel

    %%DEBUG? Why make neural event dobule wide????
    %%     U = conv(conv(Ugaussian,ones(1,2)),NEV);    % convolve neural
    %%                                                 % events with
    %%                                                 % kernel
    U = conv(Ugaussian,NEV);    % convolve neural events with kernel
    U = U(1:N)/max(U(:));                       % scale to range
                                                % [0,1]

    NEValt = U;
    
    %integrate the balloon model
    SCKS = spm_DEM_generate(M,U,{P},{6,8},{8});  %%% What are these params?
    
    % desample BOLD @ frequency observed, i.e., 1Hz ( samples per second = 1/M(1).E.dt )
    BLDgen = (SCKS.Y).';

    
    return
