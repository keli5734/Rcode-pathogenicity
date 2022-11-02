 clear 
 clc 
 close all
% 
%%

Parameter_table_HP = readtable('H1N1.HP.para.xls');
Parameter_table_LP = readtable('H1N1.LP.para.xls');
%  
% 
%  Parameter_table_HP = readtable('H5N1.HP.para.xls');
%  Parameter_table_LP = readtable('H5N1.LP.para.xls');
% 

%%

% TD_table = readtable('TD_H1N1.xls','ReadVariableNames',true);
% TD_HP = sort(TD_table{:,"HP"});
% TD_LP = sort(TD_table{:,"LP"});
% 
% TD_HP_median = TD_HP(3000);
% TD_LP_median = TD_LP(3000);
% 
% %%
% AUCD_table = readtable('AUCD_H1N1.xls','ReadVariableNames',true);
% AUCD_HP = sort(AUCD_table{:,"HP"});
% AUCD_LP = sort(AUCD_table{:,"LP"});
% 
% AUCD_HP_median = AUCD_HP(3000);
% AUCD_LP_median = AUCD_LP(3000);

%%
% load HP parameters
sV_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_1'}};
beta_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_2'}};
qFI_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_3'}};
qFM_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_4'}};
kappaA_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_5'}};
sM_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_6'}};
qprime_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_7'}};
V0_HP = 10.^Parameter_table_HP{:,{'log10_theta_HP_8'}};


% load LP parameters 
sV_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_1'}};
beta_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_2'}};
qFI_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_3'}};
qFM_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_4'}};
kappaA_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_5'}};
qprime_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_6'}};
sM_LP = sM_HP;
V0_LP = V0_HP;

len_para = length(beta_HP);%% Model parameters 
 
%%

% HP
delta_MR = 1.1e-2; % 0.011

k_m1 = 0.3; % TB paper (0.3, 2) % 0.36

k_m2 = 0.3; % TB paper

delta_M1 = 1.1e-2; % covid paper

delta_M2 = 1.1e-2;


 % Parameters in this section determine macrophage dynamics 


k1 = 0.4; %2e+3; % The maximal rate of MR -> M1 (0.4 TB paper)

k2 = 0.4*1e-5; %5e-5; % The maximal rate of MR -> M2 


V50_M = 1.2e+7; %  Half-Sat of Virus on MR --> M1 

alpha1 = 1e-4; % The effect of M2 on MR --> M1 

D_50 = 1e+6; % Half-Sat of dead cells for MR -> M2

 

% Viral infection parameters 

g_T = 0.8;

T_max = 7e+7;


%q_prime = 1e-6; %1e-5; % virus loss due to bind to macrophages

delta_I = 2;

delta_V = 5;

 

kappa_F = 3; % infected cells killed by NK cells

kappa_E =  8; %  5e-5 <--- this parameter can affect viral load



% Dead cell parameters (CoVid paper)

kappa_D = 8e-7;  % Smith et al.(2011) %8e-9; % the clearance rate of apoptotic cells by activated M1 macrophages

delta_D = 2; % degration of dead cells

 
% viral infections 

p_I =  210;

% Interferon parameters 

delta_F = 2;   % % decay rate of interferons (Ref. Cao et al. 16) 2

phi = 0.33;    % the rate of interferons turn T to R (Cao 0.33)

xi_R = 2.6;   % the rate of R turns into T (Cao 2.6)

% latent phase 

l = 4; % 4 per day which is equivalent ot 6 hours

 
% time setting of viral infection


t0_inf = 0;
t_end_inf = 10; % consider 14 days post inoculation of virus
report_point_inf = (t_end_inf - t0_inf) * 24 + 1;
report_time_inf = linspace(t0_inf, t_end_inf, report_point_inf);
len_inf = length(report_time_inf);

% Viral Variables 

MR = zeros(1,len_inf);
M1 = zeros(1,len_inf);
M2 = zeros(1,len_inf);
T = zeros(1,len_inf);
I = zeros(1,len_inf);
V = zeros(1,len_inf);
F = zeros(1,len_inf);
R = zeros(1,len_inf);
D = zeros(1,len_inf);
I_1 = zeros(1,len_inf);


options = odeset('RelTol',1e-10,'AbsTol',1e-20); 

 
%%




qFI_HP_sorted = sort(qFI_HP);
qFI_HP_median = qFI_HP_sorted(4500);
para_indx = find(qFI_HP == qFI_HP_median);

beta_base = beta_HP(para_indx);
s_M = sM_HP(para_indx);
s_V = sV_HP(para_indx);
kappa_AS = kappaA_HP(para_indx);
q_FM = qFM_HP(para_indx);
q_prime = qprime_HP(para_indx);
V_initial = V0_HP(para_indx);
q_FI = qFI_HP(para_indx);

ratio = 0:0.01:1;
beta =  (1 - ratio) * beta_base;
TD_change_HP = zeros(1, length(beta));
AUCD_change_HP = zeros(1, length(beta));


for j = 1:length(beta)

MR_inf0 = s_M / delta_MR;  % intial value of M in infection, equivalent to the number of M at day 1000 in the absence of infection
M1_inf0 = 0; % intial value of M1 in infectione, quivalent to the number of M1 at day 1000 in the absence of infection
M2_inf0 = 0;


T0 = T_max;
I0 = 0;
V0 = V_initial; 
F0 = 0;
R0 = 0;
D0 = 0;
I_10 = 0;




MR(1) = MR_inf0;
M1(1) = M1_inf0;
M2(1) = M2_inf0;
T(1) = T_max;
V(1) = V0;



init_inf = [MR_inf0,M1_inf0,M2_inf0,...
            T0,I0,V0,F0,R0,D0,I_10]';
        



[~,y_inf_HP] = ode15s(@Model, report_time_inf, init_inf, options,...
                       s_M, delta_MR,...
                       s_V, k1, k2, alpha1,D_50, k_m1, k_m2, V50_M,...
                       delta_M1,delta_M2,...
                       g_T,T_max,beta(j),phi,xi_R,delta_I,kappa_F,...
                       p_I,delta_V,...
                       q_FI,q_FM,delta_F,kappa_D,delta_D,...
                       kappa_E, kappa_AS, q_prime, l);

TD_HP(j) = (1 - min(y_inf_HP(:,4) + y_inf_HP(:,8)) / T_max) * 100; 
AUCD_HP(j) = trapz(y_inf_HP(:,9));     

TD_change_HP(j) = (TD_HP(j) - TD_HP(1)) / TD_HP(1) * 100;
AUCD_change_HP(j) = (AUCD_HP(j) - AUCD_HP(1)) / AUCD_HP(1) * 100;




end 
   

color_matrix = cell(3,1);
color_matrix{1} = [0,0,0];
color_matrix{2} = 	[0.5,0.5,0.5];
color_matrix{3} = 	'#EDB120';

LineType_matrix = cell(3,1);
LineType_matrix{1} = '-';
LineType_matrix{2} = '--';
LineType_matrix{3} = '-';

dot_plot_x = [0,0.25,0.5,0.75,1];
for i = 1:length(dot_plot_x)
dot_plot_y_indx(i) = find(ratio == dot_plot_x(i));
end 
dot_plot_y = TD_HP(dot_plot_y_indx);

figure(3) % T + R
yyaxis left
plot(ratio, TD_HP, 'Color', color_matrix{1}, 'LineWidth',3, 'LineStyle', LineType_matrix{1})
hold on 
plot(dot_plot_x,dot_plot_y,'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10)
ax = gca;
ax.FontSize = 25;
xlim([0, 1])
ax.XTick = [1,1.25,1.5,1.75,2] - 1;
xticklabels({'0','25','50','75','100'});
ax.YTick = [0,10,20,30];
yticklabels({'0','10','20','30'});
ax.FontName = 'Times new roman';
title('Max epthelium loss')
xlabel('\beta decreases (%)')
ylabel('Cell loss % (HP)')

dot_plot_y = AUCD_HP(dot_plot_y_indx);
 
 
figure(6)
yyaxis left
plot(ratio, AUCD_HP, 'Color', color_matrix{1}, 'LineWidth',3, 'LineStyle', LineType_matrix{1})
hold on 
plot(dot_plot_x,dot_plot_y,'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10)
ax = gca;
ax.FontSize = 25;
xlim([0, 1])
ax.XTick = [1,1.25,1.5,1.75,2] - 1;
xticklabels({'0','25','50','75','100'});
% ax.YTick = [-100,-80,-60,-40,-20,0];
% yticklabels({'-100%','-80%','-60%','-40%','-20%','0%'})
ax.FontName = 'Times new roman';
title('Cumulative dead cells')
xlabel('\beta decreases (%)')
ylabel('AUCD (HP)')



%%


qFI_LP_sorted = sort(qFI_LP);
qFI_LP_median = qFI_LP_sorted(4500);
para_indx = find(qFI_LP == qFI_LP_median);

beta_base = beta_LP(para_indx);
s_M = sM_LP(para_indx);
s_V = sV_LP(para_indx);
kappa_AS = kappaA_LP(para_indx);
q_FM = qFM_LP(para_indx);
q_prime = qprime_LP(para_indx);
V_initial = V0_LP(para_indx);
q_FI = qFI_LP(para_indx);


ratio = 0:0.05:1;
beta = (1 - ratio) * beta_base;
TD_change_LP = zeros(1, length(beta));
AUCD_change_LP = zeros(1, length(beta));



for j = 1:length(beta)

MR_inf0 = s_M / delta_MR;  % intial value of M in infection, equivalent to the number of M at day 1000 in the absence of infection
M1_inf0 = 0; % intial value of M1 in infectione, quivalent to the number of M1 at day 1000 in the absence of infection
M2_inf0 = 0;


T0 = T_max;
I0 = 0;
V0 = V_initial; 
F0 = 0;
R0 = 0;
D0 = 0;
I_10 = 0;




MR(1) = MR_inf0;
M1(1) = M1_inf0;
M2(1) = M2_inf0;
T(1) = T_max;
V(1) = V0;
 


init_inf = [MR_inf0,M1_inf0,M2_inf0,...
            T0,I0,V0,F0,R0,D0,I_10]';
        



[~,y_inf_LP] = ode15s(@Model, report_time_inf, init_inf, options,...
                       s_M, delta_MR,...
                       s_V, k1, k2, alpha1,D_50, k_m1, k_m2, V50_M,...
                       delta_M1,delta_M2,...
                       g_T,T_max,beta(j),phi,xi_R,delta_I,kappa_F,...
                       p_I,delta_V,...
                       q_FI,q_FM,delta_F,kappa_D,delta_D,...
                       kappa_E, kappa_AS, q_prime, l);

TD_LP(j) = (1 - min(y_inf_LP(:,4) + y_inf_LP(:,8)) / T_max) * 100; 
AUCD_LP(j) = trapz(y_inf_LP(:,9));               

TD_change_LP(j) = (TD_LP(j) - TD_LP(1)) / TD_LP(1) * 100;
AUCD_change_LP(j) = (AUCD_LP(j) - AUCD_LP(1)) / AUCD_LP(1) * 100;




end 
   

dot_plot_x = [0,0.25,0.5,0.75,1];
for i = 1:length(dot_plot_x)
dot_plot_y_indx(i) = find(ratio == dot_plot_x(i));
end 
dot_plot_y = TD_LP(dot_plot_y_indx);

figure(3) % T + R
yyaxis right
plot(ratio, TD_LP, 'Color', color_matrix{1}, 'LineWidth',3, 'LineStyle', LineType_matrix{2})
hold on 
plot(dot_plot_x,dot_plot_y,'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10)
legd = legend({'HP','','LP',''}, 'location','northeast');
legend box off
legd.FontSize = 25;
box off
ylabel('Cell loss % (LP)')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YTick = [0,0.06,0.12];
yticklabels({'0','0.06','0.12'}); 

% legend box off


dot_plot_y = AUCD_LP(dot_plot_y_indx);
 
 
figure(6)
yyaxis right 
plot(ratio, AUCD_LP, 'Color', color_matrix{1}, 'LineWidth',3, 'LineStyle', LineType_matrix{2})
hold on 
plot(dot_plot_x,dot_plot_y,'o','MarkerFaceColor','black','MarkerEdgeColor','black', 'MarkerSize',10)
box off
ylabel('AUCD (LP)')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.YTick = [0,0.5e+6,1e+6,1.5e+6,2e+6];
%yticklabels({'0','1','2','3'}); 
% legend({'HP','','LP',''}, 'location','northeast')
% legend box off




