function [K_y,G,Y_0] = SLICK(A, dt, nt_train, gamma1, gamma2, t_remove, shift, varargin )

%%  Stochastic Low-dimensional Inflated Convolutional Koopman model (SLICK)  

%       Inputs:  
%
%            A:  Convolutional expansion coeffcients obtained from tcoeffs(...)
%           dt:  Time step
%     nt_train:  Size of the training set
%    gamma 1&2:  Ridge parameters for L2 regularization
%     t_remove:  Removal of the first few snapshots; not necessary
%        shift:  Starting point of the training set

%      Outputs:  
%
%          K_y:  Inflated Koopman operator
%            G:  De-whitening filter
%          Y_0:  All the inflated state vectors

% Reference:
%        [1] T. Chu, O. T. Schmidt, Spectral reduced-order modeling of 
%        turbulent flows.      (Under preparation)
%        [2] T. Chu, O. T. Schmidt, A stochastic SPOD-Galerkin model for
%        broadband turbulent flows. Theoretical and Computational Fluid
%        Dynamics 35, no. 6 (2021): 759-782.


% T. Chu (tchu72@gatech.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision:  12-Aug-2024 Tianyi Chu <tchu72@gatech.edu>


%%


Nf    =  size(A,1);
M_n   =  size(A,2);
nt    =  size(A,3);


if nargin == 8 
    data_type  = varargin{1};
else
    data_type  = 'real';
end

if strcmpi(data_type,'real')
    
    nDFT       = (Nf-1)*2;
    
elseif strcmpi(data_type,'complex')
    
    nDFT       =  Nf;
    
end



%%   Koopman approach for convolutional coordinates, Eqns(2.24-2.25)

    X          =   reshape( permute(A(:,:,(1:nt_train-nDFT-1)+nDFT/2),[2 1 3]),[],nt_train-nDFT-1);
    Y          =   reshape( permute(A(:,:,(2:nt_train-nDFT)+nDFT/2),[2 1 3]),[],nt_train-nDFT-1);
 
    
    K_0        =   Y*(X'/(X*X'+gamma1*speye((Nf)*M_n)) ) ;
    K1         =   (K_0-eye(length(K_0)))/dt;   


 
 %% temporal derivatives
    
    B          =   zeros(M_n*(Nf),nt);
    dbdt       =   zeros(M_n*(Nf),nt);
    dadt       =   zeros(M_n*(Nf),nt);
    A_1        =  reshape(permute(A,[2 1 3]),[],size(A,3));
    
    for it=1:nt
        disp(['computing B from data at time step ' num2str(it) '/' num2str(nt)])
        
        if it<nt-1
            
            dadt(:,it)    =   (A_1(:,it+1)-A_1(:,it))/(dt);
            B(:,it)       =   (A_1(:,it+1)-A_1(:,it))/(dt)-K1*A_1(:,it);
            dbdt(:,it)    =   (A_1(:,it+2)+A_1(:,it)-2*A_1(:,it+1))/(dt^2)-K1*(A_1(:,it+1)-A_1(:,it))/(dt);
            
        elseif it==nt
            dadt(:,it)    =   (3*A_1(:,nt)+A_1(:,nt-2)-4*A_1(:,nt-1))/(2*dt);
            B(:,it)       =   (3*A_1(:,nt)+A_1(:,nt-2)-4*A_1(:,nt-1))/(2*dt)-K1*A_1(:,nt);
            dbdt(:,it)    =   (-A_1(:,it)-A_1(:,it-2)+2*A_1(:,it-1))/(dt^2)-K1*(3*A_1(:,nt)+A_1(:,nt-2)-4*A_1(:,nt-1))/(2*dt);
        end
    end
    
 %%  Inflated Koopman approach. Eqns(3.1-3.4)

 
 e_L                 =   size(X,2)-t_remove-1;         % Size of the actual training set after removal
 
 Y_0                 =   [A_1(:,(1:nt-t_remove)+t_remove); B(:,(1:nt-t_remove)+t_remove)]; 
 Y_1                 =   Y_0(:,(1:e_L)+shift);
 Y_2                 =   Y_0(:,(1:e_L)+shift+1);
    
 M                   =   dbdt(:,(1:e_L)+t_remove+shift )* (Y_1'/(Y_1*Y_1'+gamma2*speye((Nf)*M_n*2) ) );
 
 K_y                 =   [K1  eye(length(K1)); M];     % Inflated Koopman operator
 
 %% De-whitening filetr constructions, Eqns(3.5-3.12)
 
 T                   =   K_y*dt + eye(size(K_y));      % State transition matrix
    
 P_y2                =   (Y_2*Y_2'/e_L);   
 P_y1                =   (Y_1*Y_1'/e_L);
       
 R_rr                =   P_y2-T*P_y1*T';
 R_rr                =   (R_rr+R_rr')/2;
 R_rr                =   R_rr((Nf)*M_n+1:end,(Nf)*M_n+1:end );
 
[R_v,R_d]            =   eig(R_rr/dt );
R_d                  =   diag(R_d);
idx_R1               =   R_d>0;
idx_R2               =   R_d<0;
R_d_s1               =   sum(R_d);  
R_d(idx_R2)          =   0;
R_d(idx_R1)          =   R_d(idx_R1) / sum(R_d(idx_R1)) *(R_d_s1 );

G                    =   R_v(:,:)*diag(sqrt(R_d));   % de-whitening filter
    
end

