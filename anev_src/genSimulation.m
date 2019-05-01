function data = genSimulation(TS,HRFtrans, N_min, N_max, FG, LAG_min, LAG_max,A_min,A_max,B_min,B_max)

    rng('shuffle');
       
    %%Construct FN
    N_range = N_max - N_min;
    N = N_min + floor( ( N_range + 1 ) * rand( 1, 1 ) );
    
    eyeN = eye( N );
    yeyN = eyeN == 0;
   
    
    FN_min = sigmoid( A_min );
    FN_max = sigmoid( A_max );
    FN_range = FN_max - FN_min;
    FN = ( FN_min + ( FN_range * rand( N, N )) ) .* yeyN;

    FN = FN_min * yeyN;
    indices = randsample(N^2,N^2);
    [I,J] = ind2sub([N,N],indices);
    count = 0;
    i=1;
    while(count~=N)
        if(I(i)~=J(i))
            FN(I(i),J(i))=FN_max;
            count = count+1;
        end
        i=i+1;
    end
    
    %%Construct Lag Model
    LAG_range = LAG_max - LAG_min;
    LAG = (LAG_min + LAG_range * rand( N, N )) .* yeyN;
    
    %%Simumlate random external activity
    B_min = 0.05;
    B_max = 0.05;
    B_range = B_max - B_min;
    B = B_min + B_range * rand( N, 1 );
    
    %Hangle lag conversion to discrete steps
    LAGstep = floor(LAG/(1./FG));
    LAGrem = (LAG/(1./FG))-LAGstep;  %Leftover time that doesn't fit
        
    %%Simulate functional activity
    NEVgen = simulateFuncActivity( FN, LAGstep, LAGrem, TS,HRFtrans, FG, B );
    
    %%Prep simulation components for storage (clean-up)
    LAG_I = LAGstep;
    LAG_F = LAGrem;
    NEVgen = min( NEVgen, 1 );
    
    %%Save Simulation components
    data.LAG_F = LAG_F;
    data.LAG_I = LAG_I;
    data.NEVgen = NEVgen;
    data.FN = FN;
    data.N = N;
        
end