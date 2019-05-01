% Author:      Keith Bush, PhD
% Institution: University of Arkansas at Little Rock
% Date:        Aug. 12, 2013

function [encoding] = deconvolve_Bush_2013(BLDobs,kernel,nev_lr,epsilon,beta)
%-Description
% This function deconvolves the BOLD signal using Bush 2011 method
%
%-Inputs
% BLDobs - observed BOLD signal
% kernel - assumed kernel of the BOLD signal
% nev_lr - learning rate for the assignment of neural events
% epsilon - relative error change (termination condition)
% beta - slope of the sigmoid transfer function (higher = more nonlinear)
%
%-Outputs
% encoding - reconstructed neural events

  %Calc time related to observations
  N = numel(BLDobs);

  %Calc simulation steps related to simulation time
  K = numel(kernel);
  A = K-1+N;
  
  %Termination Params
  preverror = 1E9;
  currerror = 0;

  %Construct activation vector
  activation = zeros(A,1)+(2E-9).*rand(A,1)+(-1E-9);
    
  %Presolve activations to fit target_adjust as encoding
  max_hrf_id_adjust = find(kernel==max(kernel))-1;
  BLDobs_adjust = BLDobs(max_hrf_id_adjust:N);
  pre_encoding = BLDobs_adjust-min(BLDobs_adjust);
  pre_encoding = pre_encoding./max(pre_encoding);
  encoding = pre_encoding;
  activation(K:(K-1+numel(BLDobs_adjust))) = (1/beta)*log(pre_encoding./(1-pre_encoding));

  while abs(preverror-currerror) > epsilon

    %Compute encoding vector
    encoding = sigmoid(activation,beta);
    
    %Construct feature space
    feature = generate_feature(encoding,K);

    %Generate virtual bold response
    ytilde = feature(K:size(feature,1),:)*kernel;

    %Convert to percent signal change
    meanCurrent = mean(ytilde);
    brf = ytilde - meanCurrent;
    brf = brf/meanCurrent;
    
    %Compute dEdbrf
    for u=size(brf,2)
        dEdbrf(:,u) = brf(:,u)-BLDobs';
    end

    %Assume normaization does not impact deriv much.
    dEdy = dEdbrf;

    %Compute dEdy
    dEdy = dEdbrf;
    dEde = kernel;
    back_error = [zeros(1,K-1),dEdy',zeros(1,K-1)]; 
      
    %Backpropagate Errors
    delta = zeros(1,A);
    for i = 1:A
      active = activation(i);
      deda = dsigmoid(active,beta);
      dEda = dEde .* deda;
      this_error = back_error(i:((i - 1) + K));
      delta(i) = sum(dEda' .* this_error);
    end

    %Update estimate
    activation = activation-nev_lr.*delta';

    %Iterate Learning
    preverror = currerror;
    currerror = sum(dEdbrf.^2);

  end

end
  
% Support functions
function y = sigmoid(x,beta)
    y=1./(1+exp(-beta*x));
end


function y = dsigmoid(x,beta)
    y=(1-sigmoid(x,beta)).*sigmoid(x,beta); 
end

function fmatrix = generate_feature(encoding,K)
    encoding = reshape(encoding,numel(encoding),1);
    fmatrix = zeros(numel(encoding),K);
    fmatrix(:,1) = encoding;

    for i = 2:K
        fmatrix(:,i) = [zeros(i-1,1);encoding(1:(numel(encoding)-(i-1)),1)];
    end
end

