% step 1: run rbc_log.mod
% step 2: run nk_basic.mod
% step 3: run nk_wage_sticky.mod
% step 4: run this script

variables = {'Y', 'I', 'C', 'R', 'K', 'W', 'L', 'A'};
horizon = 1:40; % Adjust if necessary

figure;
for i = 1:length(variables)
    subplot(3,3,i);
    plot(horizon, rbc_results.([variables{i} '_e'])-rbc_results.P_e, 'b', 'LineWidth', 1); % RBC in blue
    hold on;
    plot(horizon, nk_results.([variables{i} '_e'])-nk_results.P_e, 'r', 'LineWidth', 1); % NK in red
    hold on;
    plot(horizon, nk_sticky_results.([variables{i} '_e'])-nk_sticky_results.P_e, 'Color', [0, 0.5, 0], 'LineWidth', 1); %Nk with sticky
    title(variables{i});
    yline(0, 'k', 'LineWidth', 1); % Zero line (dotted)
    grid on;
end
legend('RBC', 'NK', 'NK sticky');
