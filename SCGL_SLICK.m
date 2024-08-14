% SLICK model for the stochastic complex Ginzburg-Landau equation  


% Reference:
%        [1] T. Chu, O. T. Schmidt, Stochastic reduced-order Koopman model 
%        for turbulent flows.      (Under preparation)


% T. Chu (tchu72@gatech.edu; tic173@ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
% Last revision:  12-Aug-2024 Tianyi Chu <tchu72@gatech.edu; tic173@ucsd.edu>



%%  Read data

% import data

load('data_SCGL_nonlinear_correlatedforcing')

q_m         =      data{1};
dt          =      data{2};
X           =      data{3};


%%  SPOD

nDFT        =       32; 
Nf          =       nDFT;
novlp       =       Nf*3/4;
nt_train    =       20000;
nt          =       size(q_m,2);
N           =       length(X);
window      =       hamming(nDFT); 
winWeight   =       1/mean(window);

[L,P,f1]    =       spod( permute( (q_m(:,(1:nt_train))) ,[2 1]),window,[],novlp,dt);  
 
%%   convolutional expansion coefficients

M_n          =      2;  % model rank
P1           =      reshape(squeeze(P),Nf,N,[]);
[A]          =      tcoeffs( (transpose(q_m(:,1:nt))),P1,window,[],M_n)*nDFT/winWeight;  % convolutional expansion coefficients
[A1]         =      invtcoeffs(A,window,'firsthalf','complex');                          % pre-weighted convolutional expansion coefficients

if M_n>1   
    Psi_1    =      reshape(permute(reshape(squeeze(P(:,:,1:M_n)),(Nf),N,M_n),[2  3  1]),N,[]);  % reshape the SPOD modes
elseif M_n==1
    Psi_1    =      reshape(permute(reshape(squeeze(P(:,:,1:M_n)),(Nf),N),[2  1]),N,[]);
end    

q_c          =      (Psi_1*A1);                   %  reconstruction of the flow field

plot_overview;

 %%    The stochastic low-dimensional inflated convolutional Koopman (SLICK) model  


 gamma1              =   0;                         %   Ridge parameter #1
 gamma2              =   0;                         %   Ridge parameter #2
 t_remove            =   25;                        %   Removal of the dataset; not necessary
 shift               =   100;                       %   Starting point of the training set
 
[K_y,G,Y_0]          =   SLICK(A, dt, nt_train, gamma1, gamma2, t_remove, shift, 'complex' );

 
 %%  Simulations for hindcast
 
 k_max               =   500;                       %   Total number of Monte-Carlo simulations
 it0                 =   500+nDFT/2;                %   Initial condition: a random data point within the training set
 y0                  =   Y_0(:,it0+shift);
 t_span              =   400;
 
 [Y_mc]              =   SLICK_simulation(K_y, G, y0, dt, Nf, M_n,  t_span, k_max);
 [A_data]            =   invtcoeffs( permute(reshape((Y_0(1:(Nf)*M_n,(1: t_span)+it0-1+shift)),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');

 clear A_model;

for k = 1:k_max
     A_model(:,:,k)  =   invtcoeffs( permute(reshape( squeeze(Y_mc(k,1:(Nf)*M_n,(1: t_span))),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');
end

 plot_SLICK;

 SPOD_check          =   false;    % check the SPOD eigenvalues for long-time statistics: t_span = nt.

 %%  Simulations for forecast
 
 k_max               =   200;                        %   Total number of Monte-Carlo simulations
 it0                 =   size(X,2)-t_remove;         %   Initial condition: last data point of the test set
 y0                  =   Y_0(:,it0+shift);
 t_span              =   400;
 [Y_mc]              =   SLICK_simulation(K_y, G, y0, dt, Nf, M_n, t_span, k_max);
 [A_data]            =   invtcoeffs( permute(reshape((Y_0(1:(Nf)*M_n,(1:t_span)+it0-1+shift)),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');

 clear A_model;
 
for k = 1:k_max
     A_model(:,:,k)  =   invtcoeffs( permute(reshape( squeeze(Y_mc(k,1:(Nf)*M_n,(1:t_span))),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');
end

 SPOD_check          =   false;
 plot_SLICK; 
    