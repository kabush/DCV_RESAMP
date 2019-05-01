function [error] = compute_BLD_error(variant,signal)

%%    for(i=1:size(variant,1))
%%
%%        %
%%        variant_adj(i,:) = variant(i,:)-min(variant(i,:));
%%        variant_adj(i,:) = variant_adj(i,:)/max(variant_adj(i,:));
%%
%%        %
%%        signal_adj(i,:) = signal(i,:)-min(signal(i,:));
%%        signal_adj(i,:) = signal_adj(i,:)/max(signal_adj(i,:));
%%    
%%    end
    
    error = variant-signal;
    error = sqrt(mean(error.^2,2));

end



    
