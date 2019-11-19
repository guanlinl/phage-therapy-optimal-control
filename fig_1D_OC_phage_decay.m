% read data
minimal_dosage_data_1D_OC = ...
    importdata('minimal_dosage_data_1D_OC_phage_decay.txt');

% phage decay range read
phage_decay_range = minimal_dosage_data_1D_OC(:,1);

% minimal dosage data read
minimal_dosage_1D_OC = minimal_dosage_data_1D_OC(:,2);

% plot
figure(1); 
loglog(phage_decay_range(minimal_dosage_1D_OC > 0),minimal_dosage_1D_OC(minimal_dosage_1D_OC > 0),...
    '-ok','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g'); hold on;
loglog(phage_decay_range(1),minimal_dosage_1D_OC(1),...
    '+k','MarkerSize',15,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r'); hold on;
loglog(phage_decay_range(13),minimal_dosage_1D_OC(13),...
    '*k','MarkerSize',15,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b');

legend('1D-OC therapy', '1D-OC example (slow phage decay)',...
    '1D-OC example (fast phage decay)');
h = legend('1D-OC therapy', '1D-OC example (slow phage decay)',...
    '1D-OC example (fast phage decay)');
set(h,'FontName','Times New Roman','FontSize',16, 'Interpreter','latex');
legend boxoff;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex'); % set Xtick/Ytick in latex interp
yticks([1e2 1e4 1e6 1e8 1e10]); 
xticks([1e-2 1e-1 1 1e1]); 
xlabel('Phage decay rate, $$\omega\ (h^{-1})$$', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex'); 
ylabel('Dosage (PFU/g)', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex');
axis square;


% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 15; Height = 15;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

Name = '1D-phage_decay';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');


