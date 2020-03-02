% read data
treatment_data_2D_OC = ...
    importdata('minimal_dosage_data_2D_OC_tf.txt');
treatment_data_2D_OC_practical = ...
    importdata('minimal_dosage_data_2D_OC_tf_practical.txt');
% tf_range read
tf_set = treatment_data_2D_OC.data(:,1);

% minimal dosage data read
dosage_PS = treatment_data_2D_OC.data(:,2);
% minimal dosage data read
dosage_PR = treatment_data_2D_OC.data(:,3);
% minimal dosage data read
timing_PS = treatment_data_2D_OC.data(:,4);
% minimal dosage data rea
timing_PR = treatment_data_2D_OC.data(:,5);

dosage_PS = treatment_data_2D_OC_practical.data(:,2);
dosage_PR = treatment_data_2D_OC_practical.data(:,3);

% plot
figure(1); 
subplot(1,2,1);
semilogy(tf_set, dosage_PS,...
    '-or','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r'); hold on
semilogy(tf_set, dosage_PR ,...
    '->k','MarkerSize',13,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b'); 

legend('dosage of phage $P_{S}$', 'dosage of phage $P_{R}$');
h = legend('dosage of phage $P_{S}$', 'dosage of phage $P_{R}$');
set(h,'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
legend boxoff;
set(gca,'FontSize',20); set(gca,'TickLabelInterpreter', 'latex');
axis([min(tf_set), max(tf_set), 1e2, 1e10]);
yticks([1e2 1e4 1e6 1e8 1e10]);
xticks(min(tf_set): 12: max(tf_set)); 
xlabel('Final time, $$t_{f}$$ ', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex'); 
ylabel('Dosage (PFU/g)', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex');
axis square; %axis tight;

subplot(1,2,2);
plot(tf_set,timing_PS,...
    'or','MarkerSize',15); hold on
plot(tf_set,timing_PR,...
    '->b','MarkerSize',13); 
legend('injection timing of phage $P_{S}$', 'injection timing of phage $P_{R}$');
h = legend('injection timing of phage $P_{S}$', 'injection timing of phage $P_{R}$');
set(h,'FontName', 'Times New Roman','FontSize',15, 'Interpreter','latex');
legend boxoff;
set(gca,'FontSize',20);set(gca,'TickLabelInterpreter', 'latex');
axis([min(tf_set), max(tf_set), 0, 10]);
yticks(0: 2: 10);
xticks(min(tf_set): 12: max(tf_set)); 
xlabel('Final time, $$t_{f}$$ ', 'FontName', 'Times New Roman','FontSize',20,'Interpreter','latex'); 
ylabel('Hours post infection', 'FontName', 'Times New Roman','FontSize',20, 'Interpreter','latex');
axis square; %axis tight;

% output figures in eps 
ff = figure(1);
ff.Units = 'inches';
Width = 30; Height = 15;

ff.PaperSize = [Width, Height];
ff.PaperPosition = [0 0 Width, Height];
ff.Position = [0 0 Width, Height];

Name = '2d_FinalTime';
print(ff, Name, '-depsc2','-r600');
print(ff, Name, '-dpdf','-r600');