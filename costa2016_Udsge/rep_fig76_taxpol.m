% step 1: run nk_gov_fiscal_adj_1.mod to nk_gov_fiscal_adj_4.mod
% step 2: run this script

variables = {'Y', 'L', 'C', 'G', 'IG', 'TRANS', 'T', 'B'};
horizon = 1:40; % Adjust if necessary

figure;
for i = 1:length(variables)
    subplot(3,3,i);
    plot(horizon, nk_gov_res_adj1.([variables{i} '_e']), 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1); %0:black
    hold on;
    plot(horizon, nk_gov_res_adj2.([variables{i} '_e']), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1); % 0.1:blue
    hold on;
    plot(horizon, nk_gov_res_adj3.([variables{i} '_e']), 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1); %1:red
    hold on;
    plot(horizon, nk_gov_res_adj4.([variables{i} '_e']), 'Color', [0, 0.5, 0], 'LineWidth', 1); %10:green
    title(variables{i});
    yline(0, 'k', 'LineWidth', 1); % Zero line
    grid on;
end
legend('0', '0.1', '1','10');
