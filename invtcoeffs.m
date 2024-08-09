function [A1] = invtcoeffs(A,window,varargin)
%INVTCOEFFS  Inversion (or pre-weighting) of spod-based continuously-discrete temporal expansion coefficients


%   [A1] = INVTCOEFFS(A,WINDOW,VARARGIN) returns the
%   inversion of continuously-discrete temporal SPOD mode expansion coefficients of the
%   leading NMODES modes, and is performed in the low-dimensional space.
%   A is continuously-discrete temporal expansion coefficients of SPOD modes 
%   returned by [A] = tcoeffs(X,P,window,weight,nModes).
%   WINDOW is the same variable as for SPOD. 

%   reconst_opt: perform the inversion based on the center expansion
%   coefficet/full window/first half of the window  ['center'|'full'|'first half']
%   If not specified, reconst_opt defaults to 'center'.

%   data_type: complex-valuedity of X ['real'|'complex']
%   If not specified, X defaults to be real-valued.


%   The full data matrix X can be reconstructed as
%   X = reshape(permute(reshape(squeeze(P(:,:,1:nModes)),size(A,1),size(P,2),size(A,2)),[2 3 1]),size(P,2),[])*A1;
%
%   Reference:
%        [1] T. Chu, O. T. Schmidt, A stochastic SPOD-Koopman two-level
%        model for turbulent flows.      (Under preparation)


% T. Chu (tchu72@gatech.edu; tic173@ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision:  04-Aug-2023 Tianyi Chu <tic173@ucsd.edu>



dims        = size(A);
nt          = dims(3);
nFreq       = dims(1);
nModes      = dims(2);


if nargin==3
    reconst_opt = varargin{1};
else
    reconst_opt = 'center';
end

if nargin==4
    data_type = varargin{2};
else
    data_type = 'real';
end



window = window(:);



% default window size and type
if length(window)==1
    window  = hammwin(window);
end
nDFT        = length(window);
winWeight   = 1/mean(window);


A           = reshape(permute(A,[2 1 3]),[],nt);


winCorr_fac= winWeight/nDFT;


disp(' ')
disp('Calculating expansion coefficients')
disp('------------------------------------')

A1 = [];



for it = 1:nt
    
    disp(['time ' num2str(it) '/' num2str(nt)])
    

    
    if strcmpi(reconst_opt,'center')
        
        if strcmpi(data_type,'real')
            
            A1(:,it)  = transpose([ones(1,nModes) reshape([2*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+0)*([1:nDFT/2]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT/2+1)*nModes, it))/nDFT;
            
        elseif strcmpi(data_type,'complex')
            
            A1(:,it)  = transpose([ reshape([1*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+0)*([0:nDFT-1]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT)*nModes, it))/nDFT;
            
        end
        
    end
    
    
    if strcmpi(reconst_opt,'full')
        
        if strcmpi(data_type,'real')
            
            if it >= nDFT/2  &&  it <nt-nDFT/2+1
                
                Aj = 0;
                
                for DT = -nDFT/2:nDFT/2-1
                    
                    Aj           =     Aj +transpose([ones(1,nModes) reshape([2*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+DT)*([1:nDFT/2]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT/2+1)*nModes, it-DT))/nDFT*winCorr_fac;
                    
                end
                
                A1(:,it)  = Aj;
                
            else
                
                A1(:,it)         =      transpose([ones(1,nModes) reshape([2*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+0)*([1:nDFT/2]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT/2+1)*nModes, it))/nDFT;
                
            end
            
        elseif strcmpi(data_type,'complex')
            
            if it >= nDFT/2  &&  it <nt-nDFT/2+1
                
                Aj = 0;
                
                for DT = -nDFT/2:nDFT/2-1
                    
                    Aj           =     Aj +transpose([ reshape([1*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+DT)*([0:nDFT-1]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT)*nModes, it-DT))/nDFT*winCorr_fac;
                    
                end
                
                A1(:,it)  = Aj;
                
            else
                
                A1(:,it)         =      transpose([ reshape([1*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+0)*([0:nDFT-1]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT)*nModes, it))/nDFT;
                
            end
        end
        
    end
    
    
    if strcmpi(reconst_opt,'firsthalf')
        
        if strcmpi(data_type,'real')
            
            Aj = 0;
            
            for DT = 0:min(nDFT/2-1,it-1)
                
                Aj           =     Aj +transpose([ones(1,nModes) reshape([2*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+DT)*([1:nDFT/2]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT/2+1)*nModes, it-DT))/nDFT*winCorr_fac;
                
            end
            
            A1(:,it)  = Aj*sum(window)/sum(window(nDFT/2+2-(1:length(0:min(nDFT/2-1,it-1)))));
            
        elseif strcmpi(data_type,'complex')
            
            Aj = 0;
            
            for DT = 0:min(nDFT/2-1,it-1)
                
                Aj           =     Aj +transpose([ reshape([1*ones(nModes,1)]* exp(1i*2*pi*(nDFT/2+DT)*([0:nDFT-1]/nDFT)),1,[] )]).*squeeze(A(1:(nDFT)*nModes, it-DT))/nDFT*winCorr_fac;
                
            end
            
            
            A1(:,it)  = Aj*sum(window)/sum(window(nDFT/2+2-(1:length(0:min(nDFT/2-1,it-1)))));
            
        end
        
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [window] = hammwin(N)
        %HAMMWIN Standard Hamming window of lenght N
        window = 0.54-0.46*cos(2*pi*(0:N-1)/(N-1))';
    end