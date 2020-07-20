function plot_impedance(freq, impedance, config, trast, color)

subplot(211)

semilogx(freq, real(impedance),trast,'color',color, 'linew', 2); hold on;
grid on;
ylabel('Parte real ($\frac{Pa s}{m}$)','Interpreter','latex')
set(gca, 'xtick', config.plotx); set(gca,'xticklabel', config.strx); 
ylim([-20 5]); xlim([config.fmin config.fmax])
xlabel('Frequência (Hz)');
set(gca,'fontsize', 14)
ylim([-5 20])
subplot(212)

semilogx(freq, imag(impedance),trast,'color',color, 'linew', 2); hold on;
grid on

set(gca, 'xtick', config.plotx); set(gca,'xticklabel', config.strx); 
ylim([-20 5]); xlim([config.fmin config.fmax])
xlabel('Frequência (Hz)');
ylabel('Parte imaginaria ($\frac{Pa s}{m})$','Interpreter','latex')

set(gca,'fontsize', 14)