

norm_q = [];

for it = 1:size(q_m,2)
    
    norm_q(it) = norm(q_m(:,it));
    
end


norm_q_c = [];

for it = 1:size(q_c,2)
    
    norm_q_c(it) = norm(q_c(:,it));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1);



subplot(3,4, 12)

semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L(nDFT/2+1:end,1:3); L(1:nDFT/2,1:3)],'k-')


xlabel('$\omega$')

ylabel('$\lambda$')
axis square;
axis tight

ylim([1e-2 1e2])

xlim([-3 3])

set(gca,'YAxisLocation','right')


subplot(3,4,8)


mode_idx = 1;


pcolor(ones(N,1).* ([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi), X'.* ones(1,nDFT), sqrt(mean(norm_q))*(sqrt(real([ squeeze(P(nDFT/2+1:end,:,mode_idx)); squeeze(P(1:nDFT/2,:,1)) ]).^2+imag([ squeeze(P(nDFT/2+1:end,:,mode_idx)); squeeze(P(1:nDFT/2,:,mode_idx)) ]).^2).*sqrt([L(nDFT/2+1:end,mode_idx); L(1:nDFT/2,mode_idx)]))' )

shading interp;

axis square;

axis tight

xlim([-3 3])

ylim([-30 30])

caxis([-3 3])

set(gca,'xtick',[])
set(gca,'ytick',[])

title('Energy')


ax1 = subplot(3,4,[1 2 3]);


pcolor((dt*(0:2000-1)'.*ones(1,N))', (ones(2000,1).*X)', real(q_m(:,1:2000)) ); shading interp;

caxis([-3 3])

ylabel('$x$')

set(gca,'xtick',[])


axis tight

title('data (full)')

ylim([-30 30])

xlim([0 500])


ax2 = subplot(3,4,[1 2 3]+4);


pcolor((dt*(0:2000-1)'.*ones(1,N))', (ones(2000,1).*X)', real(q_c(:,1:2000)) ); shading interp;

caxis([-3 3])

ylabel('$x$')

set(gca,'xtick',[])

title('data (low-rank)')

axis tight

ylim([-30 30])

xlim([0 500])


subplot(3,4,[1 2 3]+8)


plot(dt*(0:2000-1), real(norm_q(:,1:2000)) ,'k-'); hold on;

plot(dt*(0:2000-1), real(norm_q_c(:,1:2000)) ,'r-'); hold on;

ylabel('$\|q\|_E$')

xlabel('$t$')


axis tight

ylim([2 10])
xlim([0 500])


xlim([0 500])

legend('data (full)','data (rank 2$\times$32)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(3,4, 4)

semilogy([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi ,[L(nDFT/2+1:end,1:3); L(1:nDFT/2,1:3)],'k-')



ylabel('$\lambda$')

axis square;
axis tight
title('SPOD eigens')

ylim([1e-2 1e2])

xlim([-3 3])
set(gca,'xtick',[])
set(gca,'YAxisLocation','right')


ax3 = subplot(3,4,8);



mode_idx = 1;


pcolor(ones(N,1).* ([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi), X'.* ones(1,nDFT), sqrt(mean(norm_q))*(sqrt(real([ squeeze(P(nDFT/2+1:end,:,mode_idx)); squeeze(P(1:nDFT/2,:,1)) ]).^2+imag([ squeeze(P(nDFT/2+1:end,:,mode_idx)); squeeze(P(1:nDFT/2,:,mode_idx)) ]).^2).*sqrt([L(nDFT/2+1:end,mode_idx); L(1:nDFT/2,mode_idx)]))' )

shading interp;

axis square;

axis tight

xlim([-3 3])

ylim([-30 30])

caxis([-3 3])

title('Leading modes')

set(gca,'xtick',[])
set(gca,'ytick',[])

ax4 =  subplot(3,4,12);



pcolor(ones(size(L,2),1).* ([f1(nDFT/2+1:end) f1(1:nDFT/2)]*2*pi), ((1:size(L,2))-0.5)'.* ones(1,nDFT), ([L(nDFT/2+1:end,:); L(1:nDFT/2,:)]./sum([L(nDFT/2+1:end,:); L(1:nDFT/2,:)] ,2))'  );


ylabel('mode')
xlabel('$\omega$')
title('Modal energy')

ylim([0.5 6.5])

yticks([1 2 3 4 5 6])

xlim([-3 3])

caxis([-.5 .5])

set(gca,'YAxisLocation','right')


%   cbar(2) = colorbar;
%
% set(cbar(2),'YLim',[0 0.5]);
%
% cbar(2).Label.String = '\lambda^{(i)}/\Sigma_{i}\lambda^{(i)}';


  