function pp = converge_control_figure(N_iter, J_saver,Theta_saver)
figure;
index_plt = find(J_saver > 0);
subplot(1,2,1);
plot(index_plt,J_saver(index_plt),'.-k','MarkerSize',20);
xlabel('$k$','Interpreter','latex','FontSize',20);
ylabel('$\mathcal{J}(u_{k})$','Interpreter','latex','FontSize',25);
axis square;
subplot(1,2,2);
semilogy(index_plt,abs(Theta_saver(index_plt)),'.-r','MarkerSize',20);
xlabel('$k$','Interpreter','latex','FontSize',20);
ylabel('$|\Theta(u_{k})|$','Interpreter','latex','FontSize',25);
axis square;
end