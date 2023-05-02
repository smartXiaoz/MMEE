clear all;
close all;

% 打开 fig 文件
% fig = openfig('5.fig', 'invisible'); 7b
% fig = openfig('6.fig', 'invisible'); 7c
fig = openfig('8.fig', 'invisible'); 
lh= findall(gca,'type','line')
% 找到图形对象，根据需求自动调整参数
lines = findobj(fig, 'type', 'line');

% 获取数据
xdata = get(lines, 'XData');
ydata = get(lines, 'YData');

err_mee1 = ydata{8}(1,1:986);
err_mee2 = ydata{7}(1,1:986);
err_mee3 = ydata{6}(1,1:986);
err_mee4 = ydata{5}(1,1:986);
err_hampel1 = ydata{4}(1,1:986);
err_hampel2 = ydata{3}(1,1:986);
err_hampel3 = ydata{2}(1,1:986);
err_hampel4 = ydata{1}(1,1:986);


close(fig)
LEN = 1000;
figure(1),hold on;
box on;
wid = 3;
MarkerSize = 8;
Indices=50;
plot(err_mee1, '-c', 'LineWidth', wid)
plot(err_mee2, '-r', 'LineWidth', wid)
plot(err_mee3, '-g', 'LineWidth', wid)
plot(err_mee4, '-b', 'LineWidth', wid)
plot(err_hampel1, '-c', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'c','MarkerEdgeColor','c', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel2, '-r', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'r','MarkerEdgeColor','r', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel3, '-g', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'g','MarkerEdgeColor','g', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel4, '-b', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b','MarkerEdgeColor','b', 'MarkerIndices', 1:Indices:LEN)

%% LL
% h1 = '(MEE1 (\eta=0.6, {\bfL=5}, {\sigma=2}))';
% h2 = '(MEE2 ({\eta=0.68}, {\bfL=50}, {\sigma=2}))';
% h3 = '(MEE3 ({\eta=0.87}, {\bfL=100}, {\sigma=2}))';
% h4 = '(MEE4 ({\eta=0.88}, {\bfL=200}, {\sigma=2}))';
% h5 = '(MMEE-Hample1 ({\eta=0.72},  {\bfL=5}, {\sigma=2}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h6 = '(MMEE-Hample2 ({\eta=0.732}, {\bfL=10}, {\sigma=2}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h7 = '(MMEE-Hample3 ({\eta=0.954}, {\bfL=25}, {\sigma=2}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h8 = '(MMEE-Hample4 ({\eta=0.832}, {\bfL=80}, {\sigma=2}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% set(gca,'FontName','Times New Roman','FontSize',24)
% 
% h8 = '(MMEE-Hample4 ({\eta=0.3}, L=15, {\bf\sigma=0.8}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h7 = '(MMEE-Hample3 ({\eta=0.4}, L=15, {\bf\sigma=0.6}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h6 = '(MMEE-Hample2 ({\eta=0.6}, L=15, {\bf\sigma=0.4}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h5 = '(MMEE-Hample1 ({\eta=0.6},  L=15, {\bf\sigma=0.2}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h4 = '(MEE4 ({\eta=1.1}, L=15, {\bf\sigma=0.8}))';
% h3 = '(MEE3 ({\eta=1.5}, L=15, {\bf\sigma=0.6}))';
% h2 = '(MEE2 ({\eta=1.8}, L=15, {\bf\sigma=0.4}))';
% h1 = '(MEE1 (\eta=2, L=15, {\bf\sigma=0.2}))';
% set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal')


% h8 = '(MMEE-Hample4 ({\eta=0.88}, L=15, {\bf\sigma=3}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h7 = '(MMEE-Hample3 ({\eta=0.6}, L=15, {\bf\sigma=2}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h6 = '(MMEE-Hample2 ({\eta=0.5}, L=15, {\bf\sigma=1.5}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h5 = '(MMEE-Hample1 ({\eta=0.5},  L=15, {\bf\sigma=1}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
% h4 = '(MEE4 ({\eta=0.8}, L=15, {\bf\sigma=3}))';
% h3 = '(MEE3 ({\eta=0.8}, L=15, {\bf\sigma=2}))';
% h2 = '(MEE2 ({\eta=1.5}, L=15, {\bf\sigma=1.5}))';
% h1 = '(MEE1 (\eta=2, L=15, {\bf\sigma=1}))';
% set(gca,'FontName','Times New Roman','FontSize',24,'FontWeight','normal')
% h=legend(h1,h2,h3,h4,h5,h6,h7,h8);

h8 = '(MMEE-Hample4 ({\eta=1000}, L=15, {\bf\sigma=100}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
h7 = '(MMEE-Hample3 ({\eta=35}, L=15, {\bf\sigma=20}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
h6 = '(MMEE-Hample2 ({\eta=9}, L=15, {\bf\sigma=10}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
h5 = '(MMEE-Hample1 ({\eta=2},  L=15, {\bf\sigma=5}, \Delta_1=0.5, \Delta_2=2, \Delta_3=4))';
h4 = '(MEE4 ({\eta=600}, L=15, {\bf\sigma=100}))';
h3 = '(MEE3 ({\eta=11}, L=15, {\bf\sigma=20}))';
h2 = '(MEE2 ({\eta=3}, L=15, {\bf\sigma=10}))';
h1 = '(MEE1 (\eta=1.1, L=15, {\bf\sigma=5}))';
  
set(gca,'FontName','Times New Roman','FontSize',14,'FontWeight','normal')
h=legend(h1,h2,h3,h4,h5,h6,h7,h8);
set(h,'FontName','Times New Roman','FontSize',20,'FontWeight','normal');
xlabel('Iteration');
ylabel('MSD');


