function [a] = tcoeffs(X,P,window,weight,nModes)
%TCOEFFS Continuously-discrete temporal expansion coefficients of SPOD modes


%   [A] = TCOEFFS(X,P,WINDOW,WEIGHT,NMODES) returns the
%   continuously-discrete temporal SPOD mode expansion coefficients of the
%   leading NMODES modes. P is the data matrix of SPOD modes returned by
%   SPOD. X, WINDOW and WEIGHT are the same variables as for SPOD. If no
%   windowing function is specified, a Hamming window of length WINDOW will
%   be used. If WEIGHT is empty, a uniform weighting of 1 is used.
%
%   Reference:
%     [1] A. Nekkanti, O. T. Schmidt, Frequency-time analysis, low-rank
%         reconstruction and denoising of turbulent flows using SPOD,
%         Journal of Fluid Mechanics 926, A26, 2021
%     [2] T. Chu, O. T. Schmidt, Stochastic reduced-order Koopman model 
%        for turbulent flows.      (Under preparation)
%
% A. Nekkanti (aknekkan@eng.ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Previous revision: 7-Oct-2022 Brandon Yeung <byeung@ucsd.edu>

% Last revision:     14-Aug-2023 Tianyi Chu <tchu72@gatech.edu; tic173@eng.ucsd.edu>

% Calculate the inner product between the data matrix X and the SPOD modes 
% in advance. Compute the temporal expansion coefficients for a specific 
% frequency by extracting the Fourier component of the pre-computed inner 
% product (which is low-dimensiomal) at that particular frequency.



dims        = size(X);
nt          = dims(1);
nGrid       = prod(dims(2:end));
window = window(:); weight = weight(:);


% default window size and type
if length(window)==1
    window  = hammwin(window);
end
nDFT        = length(window);
winWeight   = 1/mean(window);


% inner product weight
if isempty(weight)
    weight  = ones(nGrid,1);
end
X           = reshape(X,nt,nGrid);
X           = X-mean(X,1);              % subtract mean
ndims       = size(P);
P           = permute(P,[1 length(ndims) 2:length(ndims)-1]);

if isreal(X)
    nFreq   = ceil(nDFT/2)+1;
    
    P       = reshape(P,ceil(nDFT/2)+1,ndims(end),nGrid);
    
else
    nFreq   = nDFT;
    P       = reshape(P,nDFT,ndims(end),nGrid);
end


P1     =   reshape(permute(P(:,1:nModes,:),[3 2 1]),nGrid, nFreq*nModes);

weight     = reshape(weight,1,nGrid);


% zero-padding

X          = [zeros(ceil(nDFT/2),nGrid); X; zeros(ceil(nDFT/2),nGrid);];


% inner-product

X_proj     =      P1'* (weight'.*transpose(X));
X_proj     =      reshape(X_proj, nModes,nFreq, size(X_proj,2));


a          =      zeros(nFreq, nModes, nt);


winCorr_fac= winWeight/nDFT;


disp(' ')
disp('Calculating expansion coefficients')
disp('------------------------------------')


for i=1:nt 
    
    % correction for windowing and zero-padding
    
    if (i<ceil(nDFT/2)+1)
        corr 	= 1/(winCorr_fac*sum(window(ceil(nDFT/2)-i+1:nDFT)));
    elseif (i>nt-ceil(nDFT/2)+1)
        corr	= 1/(winCorr_fac*sum(window(1:nt+ceil(nDFT/2)-i)));
    else
        corr 	= 1;
    end
    
    % performing FFT for the inner-product coefficients
    
    for j = 1: nFreq
        
        if nModes ==1
            
            a_fft = fft( (permute(squeeze(X_proj(:,j,i:i+nDFT-1)).*window,[2 1])) );
            
            a(j,:,i)  = corr*winCorr_fac*a_fft(j) ;
            
        else
            
            a_fft = fft( (permute(squeeze(X_proj(:,j,i:i+nDFT-1)).*window',[2 1])) );
            
            a(j,:,i)  = corr*winCorr_fac*a_fft(j,:) ;
            
        end
        
        
    end  
    
    disp(['time ' num2str(i) '/' num2str(nt)])
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [window] = hammwin(N)
%HAMMWIN Standard Hamming window of lenght N
window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
end