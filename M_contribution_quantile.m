
M_HP = readtable('M_contribution_H5N1_HP.xls', 'ReadVariableNames', false);
M_LP = readtable('M_contribution_H5N1_LP.xls', 'ReadVariableNames', false);
%%
interval  = 0.05:0.05:0.95;

M_HP_quantile = quantile(table2array(M_HP), interval, 1);
M_LP_quantile = quantile(table2array(M_LP), interval, 1);

%%
writematrix(M_HP_quantile','M_HP_quantile_H5N1.xls') 
writematrix(M_LP_quantile','M_LP_quantile_H5N1.xls') 


%%


phi_HP = readtable('ratio_phi_H1N1_HP.xls', 'ReadVariableNames', false);
phi_LP = readtable('ratio_phi_H1N1_LP.xls', 'ReadVariableNames', false);

kappaF_HP = readtable('ratio_kappaF_H1N1_HP.xls', 'ReadVariableNames', false);
kappaF_LP = readtable('ratio_kappaF_H1N1_LP.xls', 'ReadVariableNames', false);
%%
interval  = 0.05:0.05:0.95;

phi_HP_quantile = quantile(table2array(phi_HP), interval, 1);
phi_LP_quantile = quantile(table2array(phi_LP), interval, 1);


kappaF_HP_quantile = quantile(table2array(kappaF_HP), interval, 1);
kappaF_LP_quantile = quantile(table2array(kappaF_LP), interval, 1);
%%
writematrix(phi_HP_quantile','phi_HP_quantile.xls') 
writematrix(phi_LP_quantile','phi_LP_quantile.xls') 

writematrix(kappaF_HP_quantile','kappaF_HP_quantile.xls') 
writematrix(kappaF_LP_quantile','kappaF_LP_quantile.xls') 


%%


Rt_HP = readtable('Rt_H1N1_HP.xls', 'ReadVariableNames', false);
Rt_LP = readtable('Rt_H1N1_LP.xls', 'ReadVariableNames', false);


interval  = 0.05:0.05:0.95;

Rt_HP_quantile = quantile(table2array(Rt_HP), interval, 1);
Rt_LP_quantile = quantile(table2array(Rt_LP), interval, 1);

Rt_HP_quantile = Rt_HP_quantile .^(1/3);
Rt_LP_quantile = Rt_LP_quantile .^(1/3);

writematrix(Rt_HP_quantile','Rt_HP_quantile.xls') 
writematrix(Rt_LP_quantile','Rt_LP_quantile.xls') 





%%


TDt_HP = readtable('TDt_H1N1_HP.xls', 'ReadVariableNames', false);
TDt_LP = readtable('TDt_H1N1_LP.xls', 'ReadVariableNames', false);


interval  = 0.05:0.05:0.95;

TDt_HP_quantile = quantile(table2array(TDt_HP), interval, 1);
TDt_LP_quantile = quantile(table2array(TDt_LP), interval, 1);


writematrix(TDt_HP_quantile','TDt_HP_quantile.xls') 
writematrix(TDt_LP_quantile','TDt_LP_quantile.xls') 
