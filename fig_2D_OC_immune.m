% read data
minimal_dosage_data_2D_OC = ...
    importdata('minimal_dosage_data_2D_OC_immnue_level.txt');
pratical_data_SDs = ...
    importdata('pratical_dosage_2D_SDs_immnue_level.txt');
pratical_data_SDs_detail = ...
    importdata('data_summary_heuristic_SDs_immune_level.txt');

% immune range read
immune_range_2d = minimal_dosage_data_2D_OC(:,1);

% minimal dosage data read
minimal_dosage_2D_OC = minimal_dosage_data_2D_OC(:,2);
pratical_dosage_SDs = (1 + pratical_data_SDs(:,2)).*minimal_dosage_2D_OC;

% dosage of two types of phage
dosage_PS = pratical_data_SDs_detail(:,4);
dosage_PR = pratical_data_SDs_detail(:,5);

% plot
figure(1); 
subplot(1,2,1);
semilogy(immune_range_2d,minimal_dosage_2D_OC,...
    '-or','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r'); hold on
semilogy(immune_range_2d,pratical_dosage_SDs ,...
    '->k','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g'); 
loglog(immune_range_2d(1),minimal_dosage_2D_OC(1),...
    '+k','MarkerSize',15,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r'); hold on;
loglog(immune_range_2d(end-1),minimal_dosage_2D_OC(end-1),...
    '*k','MarkerSize',15,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b');

legend('2D-OC therapy', 'practical therapy guided by 2D-OC', ...
    '2D-OC example (low immune density)',...
    '2D-OC example (high immune density)','Location','southwest');
h = legend('2D-OC therapy', 'practical therapy (guided by 2D-OC)', ...
    '2D-OC example (low immune density)',...
    '2D-OC example (high immune density)','Location','southwest');
set(h,'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
legend boxoff;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
axis([3e6, 9e6, 1e2, 1e10]);
yticks([1e2 1e4 1e6 1e8 1e10]);
xticks(3e6:1e6:9e6); 
xlabel('Immune density, $$I_{0}$$  (cell/g)', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex'); 
ylabel('Dosage (PFU/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex');
axis square; %axis tight;

subplot(1,2,2);
semilogy(immune_range_2d,dosage_PS,...
    'dr','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r'); hold on
semilogy(immune_range_2d,dosage_PR ,...
    'dk','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b'); 
legend('dosage of phage $P_{S}$', 'dosage of phage $P_{R}$');
h = legend('dosage of phage $P_{S}$', 'dosage of phage $P_{R}$');
set(h,'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
legend boxoff;
set(gca,'FontSize',20);set(gca,'TickLabelInterpreter', 'latex');
axis([3e6, 9e6, 1e2, 1e10]);
yticks([1e2 1e4 1e6 1e8 1e10]);
xticks(3e6:1e6:9e6); 
xlabel('Immune density, $$I_{0}$$  (cell/g)', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex'); 
ylabel('Dosage (PFU/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex');
axis square; %axis tight;

% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 30; Height = 15;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

Name = '2d_immune';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');
