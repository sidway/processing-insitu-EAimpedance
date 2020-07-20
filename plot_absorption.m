%% Plot absorção
function plot_absorption(freq, absorption, config, trast, color)


semilogx(freq, absorption,trast,'color',color, 'linew', 2); grid on

set(gca, 'xtick', config.plotx); set(gca,'xticklabel', config.strx); 
ylim([0 1]); xlim([config.fmin config.fmax])
xlabel('Frequência [Hz]');
ylabel('Coeficiente de absorção [-]')

set(gca,'fontsize', 16)