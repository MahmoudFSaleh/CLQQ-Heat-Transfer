%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name  : Mahmoud Saleh and Endre Kovacs                        %%%
% Date  : 18.03.2024                                            %%%
% E-mail: mhmodsalh84@gmail.com                                 %%%
%******************************************************************
%%%                 Plot the Result                             %%%
%%%              the Heat Conduction Equation                   %%%
%%%                 in One Dimension                            %%%
%******************************************************************
% Description:
%


%******************************************************************
% Laading the Data :
clear;
close all;
clc;
load DataOfSolution

%******************************************************************
% Plotting:

%=========================================================
% The Initial Conditions
% Reference Solution and Numerical Solution-
% - at the end of the time interval

figure (1)

plot(Node.x, U0                ,'--d', 'LineWidth', 3,'markersize',8,'markerfacecolor','w')
hold on
plot(Node.x, BenchMarklSolution,'', 'LineWidth', 12,'markersize',8,'markerfacecolor','w')
hold on
plot(Node.x, U                 ,'--o', 'LineWidth', 3,'markersize',8,'markerfacecolor','w')

title('The numerical and reference solution at the end of time interval', ...
    'FontName','TimesNewRoman','fontweight','bold','fontsize',18,'fontangle','italic')
legend(' Initial Temperature',' Reference Solution', ...
    ' CLQQ method','FontName','TimesNewRoman','fontweight','bold','fontsize',18,'fontangle','italic')

XLabel = xlabel('The Index of the Node');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');
YLabel = ylabel('Temperature');
set(YLabel, 'FontSize', 16);
set(YLabel, 'FontWeight', 'bold');
grid on

%=========================================================
% Error as a Function of time Step

figure (2)

plot(Axhstep, MaxD,'--*', 'LineWidth', 3,'markersize',10,'markerfacecolor','w')

title('Error as a Function of time Step size', ...
    'FontName','TimesNewRoman','fontweight','bold','fontsize',18,'fontangle','italic')
legend(' CNe Method','FontName','TimesNewRoman','fontweight','bold','fontsize',18,'fontangle','italic')


XLabel = xlabel('Time Step Size');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');
YLabel = ylabel('Max Error');
set(YLabel, 'FontSize', 16);
set(YLabel, 'FontWeight', 'bold');

set(gca,'xScale','log')
set(gca,'yScale','log')
grid on

%=========================================================
% Error as a function of Running time:

figure (3)

plot(RunT, MaxD,'--*', 'LineWidth', 3,'markersize',10,'markerfacecolor','w')

title('Error as a Function of Running Time', ...
    'FontName','TimesNewRoman','fontweight','bold','fontsize',18,'fontangle','italic')
legend(' CNe Method','FontName','TimesNewRoman','fontweight','bold','fontsize',18,'fontangle','italic')


XLabel = xlabel('Runnin Time');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');
YLabel = ylabel('Max Error');
set(YLabel, 'FontSize', 16);
set(YLabel, 'FontWeight', 'bold');

set(gca,'xScale','log')
set(gca,'yScale','log')
grid on