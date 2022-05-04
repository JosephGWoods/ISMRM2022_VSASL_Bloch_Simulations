function plot_VSprofile(Vcut,dv,lind,mzTag,mzCont,dmz,mzTag_laminar,mzCont_laminar,dmz_laminar)

figure('Name','Mz vs velocity','Units','normalized','Position', [0.15, 1, 0.7, 0.4]);
subplot(1,2,1); hold on;
set(gca,'FontSize',14); grid on;
plot(dv(lind), mzTag(lind), 'LineWidth', 2);
plot(dv(lind), mzCont(lind), 'LineWidth', 2);
plot(dv(lind), dmz(lind), 'LineWidth', 2);
plot([Vcut,Vcut],[-1,2],'k--', 'LineWidth', 1);
plot([-Vcut,-Vcut],[-1,2],'k--', 'LineWidth', 1);
title('Velocity response','FontSize',16);
xlabel('Velocity (cm/s)','FontSize',16);
ylabel('Mz','FontSize',16);

subplot(1,2,2); hold on;
set(gca,'FontSize',14); grid on;
plot(dv(lind), mzTag_laminar, 'LineWidth', 2);
plot(dv(lind), mzCont_laminar, 'LineWidth', 2);
plot(dv(lind), dmz_laminar, 'LineWidth', 2);
plot([Vcut,Vcut],[-1,2],'k--', 'LineWidth', 1);
plot([-Vcut,-Vcut],[-1,2],'k--', 'LineWidth', 1);
title('Laminar flow','FontSize',16);
xlabel('Mean velocity (cm/s)','FontSize',16);
ylabel('Mz','FontSize',16);
legend('Label','Control','Difference','FontSize',18,'Location','southeast');

end