function plot_transferfunction(freq, transfer_function, config, trast, color)


semilogx(freq, transfer_function,trast,'color',color, 'linew', 2); grid on

set(gca, 'xtick', config.plotx); set(gca,'xticklabel', config.strx); 
ylim([0 2]); xlim([config.fmin config.fmax])
xlabel('Frequência [Hz]');
ylabel('Função Transferencia [-]')

set(gca,'fontsize', 16)