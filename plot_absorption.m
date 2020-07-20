%% Plot absor��o
function plot_absorption(freq, absorption, config, trast, color)


semilogx(freq, absorption,trast,'color',color, 'linew', 2); grid on

set(gca, 'xtick', config.plotx); set(gca,'xticklabel', config.strx); 
ylim([0 1]); xlim([config.fmin config.fmax])
xlabel('Frequ�ncia [Hz]');
ylabel('Coeficiente de absor��o [-]')

set(gca,'fontsize', 16)