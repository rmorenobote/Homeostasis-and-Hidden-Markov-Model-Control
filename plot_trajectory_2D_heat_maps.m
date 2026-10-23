

%plots stats of the trajectories

clear;
%clf;

%colormap
n = 256;
redmap = [ones(n,1), ...
          linspace(1,0,n)', ...
          linspace(1,0,n)'];
colormap(redmap);
%colorbar;


X_HMM = load('trajectories_HMM.m');
X_var_inf = load('trajectories_var_inf.m');

X_beta_alpha = load('beta_values_alpha.m');
X_beta_prob_switch = load('beta_values_prob_switch.m');

figure(5)
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 5], ...
    'PaperUnits', 'Inches', 'PaperSize', [5, 5])

pos1 = [0.15 0.65 0.3 0.3];
subplot('Position',pos1)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
x = X_HMM(:,2);
y = X_HMM(:,3);
Nx = 10;
Ny = 10;
histogram2(x, y, ...
    0.:1:Nx+0., ...
    0.:1:Ny+0., ...
    'DisplayStyle', 'tile', ...
    'ShowEmptyBins', 'on');
axis equal tight;
%set(gca, 'YDir', 'reverse');
%xticks(1:Nx);
%yticks(1:Ny);
xticks([]);
yticks([]);
xlabel('x');
ylabel('y');
% Walls
rectangle('Position', [2 3 8 1], 'FaceColor', 'k', 'EdgeColor', 'k');
rectangle('Position', [0 3 1 1], 'FaceColor', 'k', 'EdgeColor', 'k');
%rectangle('Position', [6 5 1 1], 'FaceColor', 'y', 'EdgeColor', 'y');
%colormap redmap
cb = colorbar;
cb.Label.String = 'Number of visits';


pos2 = [0.6 0.65 0.3 0.3];
subplot('Position',pos2)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
x = X_var_inf(:,2);
y = X_var_inf(:,3);
Nx = 10;
Ny = 10;
histogram2(x, y, ...
    0.:1:Nx+0., ...
    0.:1:Ny+0., ...
    'DisplayStyle', 'tile', ...
    'ShowEmptyBins', 'on');
axis equal tight;
%set(gca, 'YDir', 'reverse');
%xticks(1:Nx);
%yticks(1:Ny);
xticks([]);
yticks([]);
xlabel('x');
ylabel('y');
% Walls
rectangle('Position', [2 3 8 1], 'FaceColor', 'k', 'EdgeColor', 'k');
rectangle('Position', [0 3 1 1], 'FaceColor', 'k', 'EdgeColor', 'k');
%rectangle('Position', [6 5 1 1], 'FaceColor', 'y', 'EdgeColor', 'y');
%colormap parula
cb = colorbar;
cb.Label.String = 'Number of visits';


pos3 = [0.15 0.2 0.3 0.3];
subplot('Position', pos3)
hold on
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman')
set(gca, 'FontName', 'Times New Roman')
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
% =========================
% HMM
% =========================
x_alpha_HMM     = X_beta_alpha(:,1).';
y_alpha_HMM     = X_beta_alpha(:,3).';
error_alpha_HMM = X_beta_alpha(:,4).';
upper_alpha_HMM = y_alpha_HMM + error_alpha_HMM;
lower_alpha_HMM = y_alpha_HMM - error_alpha_HMM;
% Shaded error region
fill([x_alpha_HMM, fliplr(x_alpha_HMM)], ...
     [upper_alpha_HMM, fliplr(lower_alpha_HMM)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);
% Mean
plot(x_alpha_HMM, y_alpha_HMM, ...
     'k-', 'LineWidth', 1);
% =========================
% VI
% =========================
x_alpha_VI     = X_beta_alpha(:,1).';
y_alpha_VI     = X_beta_alpha(:,5).';
error_alpha_VI = X_beta_alpha(:,6).';
upper_alpha_VI = y_alpha_VI + error_alpha_VI;
lower_alpha_VI = y_alpha_VI - error_alpha_VI;
% Shaded error region
fill([x_alpha_VI, fliplr(x_alpha_VI)], ...
     [upper_alpha_VI, fliplr(lower_alpha_VI)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);
% Mean
plot(x_alpha_VI, y_alpha_VI, ...
     'b-', 'LineWidth', 1);
xlabel('\alpha');
ylabel('\beta, prob desired obs');
box off
hold off



pos4 = [0.6 0.2 0.3 0.3];
subplot('Position',pos4)
hold on
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
set(gca,'fontsize',11)
set(gca, 'FontName', 'Times New Roman')
% =========================
% HMM
% =========================
x_prob_switch_HMM     = X_beta_prob_switch(:,1).';
y_prob_switch_HMM     = X_beta_prob_switch(:,3).';
error_prob_switch_HMM = X_beta_prob_switch(:,4).';
upper_prob_switch_HMM = y_prob_switch_HMM + error_prob_switch_HMM;
lower_prob_switch_HMM = y_prob_switch_HMM - error_prob_switch_HMM;
% Shaded error region
fill([x_prob_switch_HMM, fliplr(x_prob_switch_HMM)], ...
     [upper_prob_switch_HMM, fliplr(lower_prob_switch_HMM)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);
% Mean
plot(x_prob_switch_HMM, y_prob_switch_HMM, ...
     'k-', 'LineWidth', 1);
% =========================
% VI
% =========================
x_prob_switch_VI     = X_beta_prob_switch(:,1).';
y_prob_switch_VI     = X_beta_prob_switch(:,5).';
error_prob_switch_VI = X_beta_prob_switch(:,6).';
upper_prob_switch_VI = y_prob_switch_VI + error_prob_switch_VI;
lower_prob_switch_VI = y_prob_switch_VI - error_prob_switch_VI;
% Shaded error region
fill([x_prob_switch_VI, fliplr(x_prob_switch_VI)], ...
     [upper_prob_switch_VI, fliplr(lower_prob_switch_VI)], ...
     [0.8 0.8 0.8], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);
% Mean
fig_prob_switch_VI = plot(x_prob_switch_VI, y_prob_switch_VI, ...
     'b-', 'LineWidth', 1);
%legend(fig_prob_switch_VI, 'Variational Control','Box','off','FontSize', 8);
%lgd.ItemTokenSize = [5, 8];
xlabel('$p_s$', 'Interpreter', 'latex');
ylabel('\beta');
box off
hold off
% Short legend lines
annotation('line', [0.75 0.78], [0.47 0.47], ...
           'Color', 'k', 'LineWidth', 1);
annotation('line', [0.75 0.78], [0.42 0.42], ...
           'Color', 'b', 'LineWidth', 1);
% Legend text
annotation('textbox', [0.78 0.495 0.0 0.0], ...
           'String', 'HMM', ...
           'EdgeColor', 'none', ...
           'FontName', 'Times New Roman', ...
           'FontSize', 9);
annotation('textbox', [0.78 0.445 0.0 0.0], ...
           'String', 'variational', ...
           'EdgeColor', 'none', ...
           'FontName', 'Times New Roman', ...
           'FontSize', 9);



% Export as vector PDF
set(gcf, 'Renderer', 'painters');
print(gcf, 'fig1.pdf', '-dpdf', '-painters');

