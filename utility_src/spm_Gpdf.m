function f = spm_Gpdf(x,h,l)
    
%  Probability Density Function (PDF) of Gamma distribution 
%  FORMAT f = spm_Gpdf(g,h,l) 
%  
%  x - Gamma-variate   (Gamma has range [0,Inf) ) 
%  h - Shape parameter (h>0) 
%  l - Scale parameter (l>0) 
%  f - PDF of Gamma-distribution with shape & scale parameters h & l 
    
    f = exp( h*log(l) + (h-1)*log(x) + (-l*x) - log(gamma(h))); %is matlab gamma func the same?
    
end
