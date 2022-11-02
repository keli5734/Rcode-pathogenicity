clear 
clc 
close all

%%

Parameter_table_HP = readtable('H5N1.HP.para.xls');
Parameter_table_LP = readtable('H5N1.LP.para.xls');

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





sM_quantile = 0.025:0.025:0.975;
sM_HP_sorted = sort(sM_HP);
s_M = sM_HP_sorted(round(sM_quantile * len_para));
for i = 1:length(sM_quantile)
para_indx(i) = find(sM_HP == s_M(i));
end 

len_sM = length(sM_quantile);
s_V = sV_HP(para_indx);
q_prime = qprime_HP(para_indx);
beta = beta_HP(para_indx);
kappa_AS = kappaA_HP(para_indx);
q_FI = qFI_HP(para_indx);
q_FM = qFM_HP(para_indx);
V_initial = V0_HP(para_indx);
 
wb = waitbar(0,'please wait...HP');


M1_prop = zeros(len_inf,len_sM);
M2_prop = zeros(len_inf,len_sM);
MR_prop = zeros(len_inf,len_sM);

for i = 1:len_sM
    for j = 1: len_inf
MR_inf0 = s_M(i) / delta_MR;  % intial value of M in infection, equivalent to the number of M at day 1000 in the absence of infection
M1_inf0 = 0; % intial value of M1 in infectione, quivalent to the number of M1 at day 1000 in the absence of infection
M2_inf0 = 0;


T0 = T_max;
I0 = 0;
V0 = V_initial(i);
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
                       s_M(i), delta_MR,...
                       s_V(i), k1, k2, alpha1,D_50, k_m1, k_m2, V50_M,...
                       delta_M1,delta_M2,...
                       g_T,T_max,beta(i),phi,xi_R,delta_I,kappa_F,...
                       p_I,delta_V,...
                       q_FI(i),q_FM(i),delta_F,kappa_D,delta_D,...
                       kappa_E, kappa_AS(i), q_prime(i), l);
                   
            MR_prop(j,i) = y_inf_HP(j,1) / (y_inf_HP(j,1) + y_inf_HP(j,2) + y_inf_HP(j,3));
            M1_prop(j,i) = y_inf_HP(j,2) / (y_inf_HP(j,1) + y_inf_HP(j,2) + y_inf_HP(j,3));  
            M2_prop(j,i) = y_inf_HP(j,3) / (y_inf_HP(j,1) + y_inf_HP(j,2) + y_inf_HP(j,3));  
            
            
            
            
    end
    
    waitbar((i-1)/len_sM)

end 



close(wb)  

writematrix(MR_prop,'MR_prop_HP_H5N1.xls');
writematrix(M1_prop,'M1_prop_HP_H5N1.xls');
writematrix(M2_prop,'M2_prop_HP_H5N1.xls');




M1_prop = zeros(len_inf,len_sM);
M2_prop = zeros(len_inf,len_sM);
MR_prop = zeros(len_inf,len_sM);

 
sM_LP_sorted = sort(sM_LP);
s_M = sM_LP_sorted(round(sM_quantile * len_para));
for i = 1:length(sM_quantile)
para_indx(i) = find(sM_LP == s_M(i));
end 

s_V = sV_LP(para_indx);
q_prime = qprime_LP(para_indx);
beta = beta_LP(para_indx);
kappa_AS = kappaA_LP(para_indx);
q_FI = qFI_LP(para_indx);
q_FM = qFM_LP(para_indx);
V_initial = V0_LP(para_indx);
 
wb = waitbar(0,'please wait...LP');

for i = 1:len_sM
    for j = 1: len_inf
MR_inf0 = s_M(i) / delta_MR;  % intial value of M in infection, equivalent to the number of M at day 1000 in the absence of infection
M1_inf0 = 0; % intial value of M1 in infectione, quivalent to the number of M1 at day 1000 in the absence of infection
M2_inf0 = 0;


T0 = T_max;
I0 = 0;
V0 = V_initial(i);
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
                       s_M(i), delta_MR,...
                       s_V(i), k1, k2, alpha1,D_50, k_m1, k_m2, V50_M,...
                       delta_M1,delta_M2,...
                       g_T,T_max,beta(i),phi,xi_R,delta_I,kappa_F,...
                       p_I,delta_V,...
                       q_FI(i),q_FM(i),delta_F,kappa_D,delta_D,...
                       kappa_E, kappa_AS(i), q_prime(i), l);
                   
            MR_prop(j,i) = y_inf_LP(j,1) / (y_inf_LP(j,1) + y_inf_LP(j,2) + y_inf_LP(j,3));
            M1_prop(j,i) = y_inf_LP(j,2) / (y_inf_LP(j,1) + y_inf_LP(j,2) + y_inf_LP(j,3));  
            M2_prop(j,i) = y_inf_LP(j,3) / (y_inf_LP(j,1) + y_inf_LP(j,2) + y_inf_LP(j,3));  
            
            
            
            
    end
    
    waitbar((i-1)/len_sM)

end 



close(wb)  

writematrix(MR_prop,'MR_prop_LP_H5N1.xls');
writematrix(M1_prop,'M1_prop_LP_H5N1.xls');
writematrix(M2_prop,'M2_prop_LP_H5N1.xls');




