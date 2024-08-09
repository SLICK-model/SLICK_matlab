



   %% analytical UQ


mj                   =   zeros ( 2*(Nf)*M_n, 200 + 1);
mj(:,1)              =   y0;

Pj                   =   zeros(2*(Nf)*M_n, 2*(Nf)*M_n ,1200+1);
pj                   =   zeros( 2*(Nf)*M_n  ,200+1);


for it=1:200
    if mod(it,1)==0
        disp([' it= ' num2str(it) ', t=' num2str(it*0.02)])
    end
    
    mj(:,it+1)    =   (eye(length(K_y))+K_y*dt)*mj(:,it);
    
    Pj(:,:,it+1)  =   (eye(length(K_y))+K_y*dt)*squeeze(Pj(:,:,it))*(eye(length(K_y))+K_y*dt)'+ 1* [zeros((Nf)*M_n,(Nf)*2*M_n);zeros((Nf)*M_n) G*G'*dt];
    
    pj(:,it+1)    =   sqrt(diag(squeeze(Pj(:,:,it+1))));

 
end

[A_upper]        =      invtcoeffs( permute(reshape((mj(1:(Nf)*M_n,(1:200))+2*pj(1:Nf*M_n,1:200)),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');

[A_lower]        =      invtcoeffs( permute(reshape((mj(1:(Nf)*M_n,(1:200))-2*pj(1:Nf*M_n,1:200)),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');



%%

 f_j = [(32-1)*M_n+1  (4-1)*M_n+1  (29-1)*M_n+1 ];   
 
 
figure;


for j=1:3
    
    subplot(4,1,j);
    plot(((1:200)-nDFT/2)*dt,real(A_data( f_j(j),(1:200) )    ),'r-','MarkerSize',4); hold on;
    plot(((1:200)-nDFT/2)*dt,real(squeeze(A_model(f_j(j) ,(1:200)+0,1 ))),'b'); hold on;
    
    t2   = [((1:200)-nDFT/2)*dt ,fliplr(((1:200)-nDFT/2)*dt )];
    
    
    inBetween  =[ (squeeze(max(real(A_model(f_j(j),1:200,1:k_max)),[],3))), fliplr((squeeze(min(real(A_model(f_j(j) ,1:200,1:k_max)),[],3) )))];
    
    %
    fill(t2,inBetween,'b','FaceAlpha',.2,'EdgeAlpha',0)
    
    plot( ((1:200)-nDFT/2)*dt, real(A_upper(f_j(j),1:200)) ,'b--') ;
    plot( ((1:200)-nDFT/2)*dt, real(A_lower(f_j(j),1:200)) ,'b--') ;
    
    
    if j ==1
        line([0 0],[-5 5] ,'LineStyle','--','Color','k')
    else
        line([0 0],[-1 1] ,'LineStyle','--','Color','k')
    end
    
    %     legend('reduced-order data','reduced-order model')
    xlim([-5 20])
    ylabel('$\mathbf{a}_j$')
    
end


   q0  = (Psi_1*A_data(:,1:200));
  
   norm_q =  sqrt(sum(conj( q_m(:,(1:200)+t_remove+shift+it0-1) ).*q_m(:,(1:200)+t_remove+shift+it0-1),1));
   
   norm_c =  sqrt(sum(conj(q0).*q0,1));

   norm_x = [];
   
   for k = 1:k_max
       
   q1  = (Psi_1*squeeze(A_model(:,1:200,k)));
   
   norm_x(:,k) =  sqrt(sum(conj(q1).*q1,1));

   end



   
   subplot(4,1,4)
   
   plot(((1:200)-nDFT/2)*dt,real(norm_q( (1:200) )    ),'k-','MarkerSize',4); hold on;
   plot(((1:200)-nDFT/2)*dt,real(norm_c( (1:200) )    ),'r-','MarkerSize',4); hold on;
    plot(((1:200)-nDFT/2)*dt,real((norm_x((1:200)+0,2 ))),'b'); hold on;
     
    t2   = [((1:200)-nDFT/2)*dt , fliplr( ((1:200)-nDFT/2)*dt )];
    
    
    inBetween  =[ (squeeze(max(real(norm_x(1:200,:)),[],2)))', fliplr((min(real(norm_x(1:200,:)),[],2))')];
    
    
    fill(t2,inBetween,'b','FaceAlpha',.2,'EdgeAlpha',0)  
     
    line([0 0],[0 15] ,'LineStyle','--','Color','m')
    
    legend('data (full)','data (rank 2$\times$32)','model (rank 2$\times$32)')
    xlim([-5 20]) 

    ylim([2 14])
    
    xlabel('$t$')
    ylabel('$\|\mathbf{q}\|_E$')
    
    drawnow;



%%

   q0  = (Psi_1*A_data(:,1:400));

   q1  = (Psi_1*squeeze(A_model(:,1:400,2)));


 
figure;



subplot(3,1,1)


pcolor((dt*((0:400-1)-(nDFT/2))'.*ones(1,N))', (ones(400,1).*X)', real(q_m(:,(1:400)+t_remove+shift+it0-1))  ); shading interp;

caxis([-5 5])

ylabel('$x$')

set(gca,'xtick',[])

title('data (full)')

% axis square; 
axis tight

line([0*dt 0*dt],[-40 40],'Color','r','LineStyle','--')

xlim([-(nDFT/2)*dt 50])
ylim([-30 30])

subplot(3,1,2)

pcolor((dt*((0:400-1)-(nDFT/2))'.*ones(1,N))', (ones(400,1).*X)', real(q0) ); shading interp;


caxis([-5 5])

ylabel('$x$')

set(gca,'xtick',[])

title('data (reduced-order)')

% axis square; 
axis tight

line([0*dt 0*dt],[-40 40],'Color','r','LineStyle','--')

ylim([-30 30])

xlim([-(nDFT/2)*dt 50])

subplot(3,1,3)


pcolor((dt*((0:400-1)-(nDFT/2))'.*ones(1,N))', (ones(400,1).*X)',  (real(q1) ) ); shading interp;

caxis([-5 5])

ylabel('$x$')

xlabel('$t$')

title('model (reduced-order)')

% axis square; 
axis tight

line([0*dt 0*dt],[-40 40],'Color','m','LineStyle','--')

ylim([-30 30])

xlim([-(nDFT/2)*dt 50])

drawnow;



    
    
 %%    SPOD check
 
 
 if SPOD_check == true
 

[A_model_all]    =      invtcoeffs( permute(reshape( squeeze(Y_mc(1,1:(Nf)*M_n,(1:end))),M_n,Nf,[]),[2 1 3]),window,'firsthalf','complex');

q_c1             =      Psi_1*A_model_all;

Q_c1             =      permute(q_c1(:,Nf:end-Nf),[2 1]);

clear q_c1

[L_c1,~,f_c1]    =      spod((Q_c1),hamming(nDFT),[],novlp,dt);

clear Q_c1

q_c              =      Psi_1*A1;
Q_c              =      permute(q_c(:,nDFT/2:end-nDFT/2),[2 1]);

clear q_c

[L_c,~,f_c]      =      spod((Q_c),hamming(nDFT),[],novlp,dt); 
clear Q_c




 
figure;



semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L(nDFT/2+1:end,1:1); L(1:nDFT/2,1:1)],'k-'); hold on;
semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L_c(nDFT/2+1:end,1:1); L_c(1:nDFT/2,1:1)],'r-'); hold on;
semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L_c1(nDFT/2+1:end,1:1); L_c1(1:nDFT/2,1:1)],'b--'); hold on;



semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L(nDFT/2+1:end,2:M_n+1); L(1:nDFT/2,2:M_n+1)],'k-'); hold on;
semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L_c(nDFT/2+1:end,2:M_n+1); L_c(1:nDFT/2,2:M_n+1)],'r-'); hold on;
semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L_c1(nDFT/2+1:end,2:M_n+1); L_c1(1:nDFT/2,2:M_n+1)],'b--'); hold on;


legend('data (full)','data (rank $2\times32$)', 'model (rank $2\times32$)')

xlabel('$\omega$')

ylabel('$\lambda$')

xlim([-3 3])

ylim([5e-3 5e1])

 
 end
 