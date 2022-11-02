clear 
clc 
close all

%%
% Parameter_table_HP = readtable('H1N1.HP.para.xls');
% Parameter_table_LP = readtable('H1N1.LP.para.xls');
 
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


%load LP parameters 
sV_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_1'}};
beta_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_2'}};
qFI_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_3'}};
qFM_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_4'}};
kappaA_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_5'}};
qprime_LP = 10.^Parameter_table_LP{:,{'log10_theta_LP_6'}};
sM_LP = sM_HP;
V0_LP = V0_HP;

len_para = length(beta_HP);



%%

TD_table = readtable('TD_H5N1.xls','ReadVariableNames',false);
TD_HP = TD_table{:,1};
TD_LP = TD_table{:,2};


%%
AUCD_table = readtable('AUCD_H5N1.xls','ReadVariableNames',false);
AUCD_HP = AUCD_table{:,1};
AUCD_LP = AUCD_table{:,2};

%%
PRCC_var = {log10(beta_HP./beta_LP),...
            log10(qFI_HP./qFI_LP),...
            log10(qFM_HP./qFM_LP),...
            log10(sV_HP./sV_LP),...
            log10(qprime_HP./qprime_LP),...
            log10(kappaA_HP./kappaA_LP)};

for i = 1:6 % 6 parameters 

PRCC_TD(i) = my_PRCC_PLOT2(PRCC_var{i}, TD_HP./TD_LP);
PRCC_AUCD(i) = my_PRCC_PLOT2(PRCC_var{i}, AUCD_HP./AUCD_LP);


end

%%

PRCC_value = [PRCC_TD, PRCC_AUCD]';

writematrix(PRCC_value, 'PRCC_values_fitted_ratio_H5N1.xls')



%%


PRCC_var_HP = {beta_HP, qFI_HP, qFM_HP, sV_HP, qprime_HP, kappaA_HP};
PRCC_var_LP = {beta_LP, qFI_LP, qFM_LP, sV_LP, qprime_LP, kappaA_LP};



%%
PRCC_TD_HP = zeros(1,6);
PRCC_TD_LP = zeros(1,6);

PRCC_AUCD_HP = zeros(1,6);
PRCC_AUCD_LP = zeros(1,6);


for i = 1:6 % 6 parameters 

PRCC_TD_HP(i) = my_PRCC_PLOT2(PRCC_var_HP{i}, TD_HP);
PRCC_TD_LP(i) = my_PRCC_PLOT2(PRCC_var_LP{i}, TD_LP);


PRCC_AUCD_HP(i) = my_PRCC_PLOT2(PRCC_var_HP{i}, AUCD_HP);
PRCC_AUCD_LP(i) = my_PRCC_PLOT2(PRCC_var_LP{i}, AUCD_LP);
 

end 
%%

PRCC_value = [PRCC_TD_HP, PRCC_TD_LP, PRCC_AUCD_HP, PRCC_AUCD_LP]';

writematrix(PRCC_value, 'PRCC_value_H5N1.xls')



