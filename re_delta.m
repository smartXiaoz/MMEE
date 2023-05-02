clear all;
close all;
clc
% 打开 fig 文件
fig = openfig('c.fig', 'invisible');

% 找到图形对象，根据需求自动调整参数
lines = findobj(fig, 'type', 'line');

% 获取数据
xdata = get(lines, 'XData');
ydata = get(lines, 'YData');

err_mee1 = ydata{5};
err_hampel1 = ydata{4};
err_hampel2 = ydata{3};
err_hampel3 = ydata{2};
err_hampel4 = ydata{1};
% % 绘制所有曲线
% figure;
% hold on;
% for i = 1:length(lines)
%     plot(xdata{i}, ydata{i});
% end
close(fig)
LEN = 1000;
figure(1),hold on;
box on;
wid = 2.5;
MarkerSize = 6;
Indices=25;
plot(err_mee1, '-c', 'LineWidth', wid, 'Marker', 'd', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'c', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel1, '-r', 'LineWidth', wid, 'Marker', 'p', 'MarkerSize',MarkerSize , 'MarkerFaceColor', 'r', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel2, '-g', 'LineWidth', wid, 'Marker', 'h', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'g', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel3, '-b', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'b', 'MarkerIndices', 1:Indices:LEN)
plot(err_hampel4, '-m', 'LineWidth', wid, 'Marker', '*', 'MarkerSize', MarkerSize, 'MarkerFaceColor', 'm', 'MarkerIndices', 1:Indices:LEN)
%% aa
% h1 = 'MEE1 (\eta=0.45, L=15, {\sigma=2})';
% h2 = 'MMEE-Hampel1 ({\eta=0.34},  L=15, {\sigma=2}, {\bf\Delta_1=0.25}, {\Delta_2=2}, {\Delta_3=4})';
% h3 = 'MMEE-Hampel2 ({\eta=0.24},  L=15, {\sigma=2}, {\bf\Delta_1=0.5}, {\Delta_2=2}, {\Delta_3=4})';
% h4 = 'MMEE-Hampel3 ({\eta=0.29},  L=15, {\sigma=2}, {\bf\Delta_1=1}, {\Delta_2=2}, {\Delta_3=4})';
% h5 = 'MMEE-Hampel4 ({\eta=0.35},  L=15, {\sigma=2}, {\bf\Delta_1=1.5}, {\Delta_2=2}, {\Delta_3=4})';

%% bb
% h1 = 'MEE1 (\eta=0.45, L=15, {\sigma=2})';
% h2 = 'MMEE-Hampel1 ({\eta=0.34},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\bf\Delta_2=1}, {\Delta_3=4})';
% h3 = 'MMEE-Hampel2 ({\eta=0.24},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\bf\Delta_2=1.5}, {\Delta_3=4})';
% h4 = 'MMEE-Hampel3 ({\eta=0.29},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\bf\Delta_2=2}, {\Delta_3=4})';
% h5 = 'MMEE-Hampel4 ({\eta=0.35},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\bf\Delta_2=3}, {\Delta_3=4})';

%% cc
h1 = 'MEE1 (\eta=0.45, L=15, {\sigma=2})';
h2 = 'MMEE-Hampel1 ({\eta=0.34},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\Delta_2=2}, {\bf\Delta_3=4})';
h3 = 'MMEE-Hampel2 ({\eta=0.24},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\Delta_2=2}, {\bf\Delta_3=6})';
h4 = 'MMEE-Hampel3 ({\eta=0.29},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\Delta_2=2}, {\bf\Delta_3=8})';
h5 = 'MMEE-Hampel4 ({\eta=0.35},  L=15, {\sigma=2}, {\Delta_1=0.5}, {\Delta_2=2}, {\bf\Delta_3=10})';

h=legend(h1,h2,h3,h4,h5);
ylim([-35, 10]);
set(h,'FontName','Times New Roman','FontSize',15,'FontWeight','normal');
xlabel('Iteration');
ylabel('MSD');
