
% Script for plotting the results of D. suzukii life tables analysis
% Created by Luca Rossini on 10 August 2024
% Last update 30 August 2024
% e-mail: luca.rossini@ulb.be
% 
% This script is recalled by 'AnalisiStepMatlab.m'
% 
% This script plots the results of all the stages

%% Plot best fit functions - EGG development

     % For loop to evaluate the expressions

BriStore_E = [];
LogStore_E = [];
LacOneStore_E = [];
LacTwoStore_E = [];
SDMStore_E = [];
TempStore_E = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_E);

    BriStore_E = [BriStore_E, Bri_TMP];
    TempStore_E = [TempStore_E, i];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_E);
    LogStore_E = [LogStore_E, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_E);
    LacOneStore_E = [LacOneStore_E, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_E);
    LacTwoStore_E = [LacTwoStore_E, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_E);
    SDMStore_E = [SDMStore_E, SDM_TMP];

end


%% Plot best fit functions - Larva 1 development

    % For loop to evaluate the expressions

BriStore_L1 = [];
LogStore_L1 = [];
LacOneStore_L1 = [];
LacTwoStore_L1 = [];
SDMStore_L1 = [];
TempStore_L1 = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_L1);

    BriStore_L1 = [BriStore_L1, Bri_TMP];
    TempStore_L1 = [TempStore_L1, i];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_L1);
    LogStore_L1 = [LogStore_L1, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_L1);
    LacOneStore_L1 = [LacOneStore_L1, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_L1);
    LacTwoStore_L1 = [LacTwoStore_L1, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_L1);
    SDMStore_L1 = [SDMStore_L1, SDM_TMP];

end


%% Plot best fit functions - Larva 2 development

    % For loop to evaluate the expressions

BriStore_L2 = [];
LogStore_L2 = [];
LacOneStore_L2 = [];
LacTwoStore_L2 = [];
SDMStore_L2 = [];
TempStore_L2 = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_L2);

    BriStore_L2 = [BriStore_L2, Bri_TMP];
    TempStore_L2 = [TempStore_L2, i];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_L2);
    LogStore_L2 = [LogStore_L2, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_L2);
    LacOneStore_L2 = [LacOneStore_L2, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_L2);
    LacTwoStore_L2 = [LacTwoStore_L2, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_L2);
    SDMStore_L2 = [SDMStore_L2, SDM_TMP];

end


%% Plot best fit functions - Larva 3 development

    % For loop to evaluate the expressions

BriStore_L3 = [];
LogStore_L3 = [];
LacOneStore_L3 = [];
LacTwoStore_L3 = [];
SDMStore_L3 = [];
TempStore_L3 = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_L3);
    BriStore_L3 = [BriStore_L3, Bri_TMP];
    TempStore_L3 = [TempStore_L3, i];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_L3);
    LogStore_L3 = [LogStore_L3, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_L3);
    LacOneStore_L3 = [LacOneStore_L3, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_L3);
    LacTwoStore_L3 = [LacTwoStore_L3, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_L3);
    SDMStore_L3 = [SDMStore_L3, SDM_TMP];

end


%% Plot best fit functions - Pupa development

    % For loop to evaluate the expressions

BriStore_P = [];
LogStore_P = [];
LacOneStore_P = [];
LacTwoStore_P = [];
SDMStore_P = [];
TempStore_P = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_P);
    BriStore_P = [BriStore_P, Bri_TMP];
    TempStore_P = [TempStore_P, i];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_P);
    LogStore_P = [LogStore_P, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_P);
    LacOneStore_P = [LacOneStore_P, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_P);
    LacTwoStore_P = [LacTwoStore_P, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_P);
    SDMStore_P = [SDMStore_P, SDM_TMP];

end


 %% Plot best fit functions - Adult males development
 
    % For loop to evaluate the expressions
 
WangStore_AM = [];
BriStore_AM = [];
LogStore_AM = [];
LacOneStore_AM = [];
LacTwoStore_AM = [];
SDMStore_AM = [];
TempStore_AM = [];
 
for j = 0 * 100 : 4500

    i = j * 0.01;

        % Wang

    Wang_TMP = Functions.Wang_Plot(i, WangRes_AM);
    WangStore_AM = [WangStore_AM, Wang_TMP];
    TempStore_AM = [TempStore_AM, i];

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_AM);
    BriStore_AM = [BriStore_AM, Bri_TMP];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_AM);
    LogStore_AM = [LogStore_AM, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_AM);
    LacOneStore_AM = [LacOneStore_AM, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_AM);
    LacTwoStore_AM = [LacTwoStore_AM, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_AM);
    SDMStore_AM = [SDMStore_AM, SDM_TMP];
 
end
 
 
%% Plot best fit functions - Adult females development

     % For loop to evaluate the expressions

WangStore_AF = [];
BriStore_AF = [];
LogStore_AF = [];
LacOneStore_AF = [];
LacTwoStore_AF = [];
SDMStore_AF = [];
TempStore_AF = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Wang

    Wang_TMP = Functions.Wang_Plot(i, WangRes_AF);
    WangStore_AF = [WangStore_AF, Wang_TMP];
    TempStore_AF = [TempStore_AF, i];

        % Briere

    Bri_TMP = Functions.Bri_Plot(i, BriRes_AF);
    BriStore_AF = [BriStore_AF, Bri_TMP];

        % Logan

    Log_TMP = Functions.Log_Plot(i, LogRes_AF);
    LogStore_AF = [LogStore_AF, Log_TMP];

        % Lactin-1

    LacOne_TMP = Functions.LacOne_Plot(i, LacOneRes_AF);
    LacOneStore_AF = [LacOneStore_AF, LacOne_TMP];

        % Lactin-2

    LacTwo_TMP = Functions.LacTwo_Plot(i, LacTwoRes_AF);
    LacTwoStore_AF = [LacTwoStore_AF, LacTwo_TMP];

        % Sharpe and De Michele

    SDM_TMP = Functions.SDM_Plot(i, SDMRes_AF);
    SDMStore_AF = [SDMStore_AF, SDM_TMP];
 
end


%% Plot best fit functions - Wang mortalities

     % For loop to evaluate the expressions

WangStore_E = [];
WangStore_L1 = [];
WangStore_L2 = [];
WangStore_L3 = [];
TempStore_Wang = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Egg

    Wang_TMP = Functions.Wang_Plot(i, WangRes_E);
    WangStore_E = [WangStore_E, Wang_TMP];
    TempStore_Wang = [TempStore_Wang, i];

        % L1

    Wang_TMP = Functions.Wang_Plot(i, WangRes_L1);
    WangStore_L1 = [WangStore_L1, Wang_TMP];

        % L2

    Wang_TMP = Functions.Wang_Plot(i, WangRes_L2);
    WangStore_L2 = [WangStore_L2, Wang_TMP];

        % L3
    
    Wang_TMP = Functions.Wang_Plot(i, WangRes_L3);
    WangStore_L3 = [WangStore_L3, Wang_TMP];

end


%% Plot best fit functions - Exponential mortality

     % For loop to evaluate the expressions

ExpPolStore_E = [];
ExpPolStore_L1 = [];
ExpPolStore_L2 = [];
ExpPolStore_L3 = [];
TempStore_ExpPol = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Egg

    ExpPol_TMP = Functions.Exp1_Plot(i, ExpPolRes_E);
    ExpPolStore_E = [ExpPolStore_E, ExpPol_TMP];
    TempStore_ExpPol = [TempStore_ExpPol, i];

        % L1

    ExpPol_TMP = Functions.Exp1_Plot(i, ExpPolRes_L1);
    ExpPolStore_L1 = [ExpPolStore_L1, ExpPol_TMP];

        % L2

    ExpPol_TMP = Functions.Exp1_Plot(i, ExpPolRes_L2);
    ExpPolStore_L2 = [ExpPolStore_L2, ExpPol_TMP];

        % L3
    
    ExpPol_TMP = Functions.Exp1_Plot(i, ExpPolRes_L3);
    ExpPolStore_L3 = [ExpPolStore_L3, ExpPol_TMP];

end


%% Plot best fit functions - Kim and Lee mortalities

     % For loop to evaluate the expressions

KLStore_E = [];
KLStore_L1 = [];
KLStore_L2 = [];
KLStore_L3 = [];
KLStore_Wang = [];
TempStore_KL = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Egg

    KL_TMP = Functions.KL_Plot(i, KLRes_E);
    KLStore_E = [KLStore_E, KL_TMP];
    TempStore_KL = [TempStore_KL, i];

        % L1

    KL_TMP = Functions.KL_Plot(i, KLRes_L1);
    KLStore_L1 = [KLStore_L1, KL_TMP];

        % L2

    KL_TMP = Functions.KL_Plot(i, KLRes_L2);
    KLStore_L2 = [KLStore_L2, KL_TMP];

        % L3
    
    KL_TMP = Functions.KL_Plot(i, KLRes_L3);
    KLStore_L3 = [KLStore_L3, KL_TMP];

end


%% Plot best fit functions - Gaussian like fertility
     % For loop to evaluate the expressions

FertStore = [];
TempStore_Fert = [];

for j = 0 * 100 : 4500

    i = j * 0.01;

        % Egg production over temperature

    Fert_TMP = Functions.Fert_Plot(i, FertRes);
    FertStore = [FertStore, Fert_TMP];
    TempStore_Fert = [TempStore_Fert, i];

end

%% Make the plots - Fit and experimental data

    % Egg - Development

figure('Name', 'Egg')
hold on

boxchart(Temp_E, E_Rates)
plot(TempStore_E, BriStore_E, '-.', 'LineWidth', 1.25)
plot(TempStore_E, LogStore_E, '--', 'LineWidth', 1.25)
plot(TempStore_E, LacOneStore_E, ':', 'LineWidth', 1.25)
plot(TempStore_E, LacTwoStore_E, '-', 'LineWidth', 1.25)
plot(TempStore_E, SDMStore_E, 'LineWidth', 1.25)

title('Egg to L1 development')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 40]);
ylim([0 1.6]);
legend('Experimental data', 'Briere', 'Logan', ' Lactin-1', 'Lactin-2', ...
       'Sharpe and De Michele')

    % Larva 1 - Development

figure('Name', 'Larva 1')
hold on

boxchart(Temp_L1, L1_Rates)
plot(TempStore_L1, BriStore_L1, '-.', 'LineWidth', 1.25)
plot(TempStore_L1, LogStore_L1, '--', 'LineWidth', 1.25)
plot(TempStore_L1, LacOneStore_L1, ':', 'LineWidth', 1.25)
plot(TempStore_L1, LacTwoStore_L1, '-', 'LineWidth', 1.25)
plot(TempStore_L1, SDMStore_L1, 'LineWidth', 1.25)

title('L1 to L2 development')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 40]);
ylim([0 1.6]);
legend('Experimental data', 'Briere', 'Logan', ' Lactin-1', 'Lactin-2', ...
       'Sharpe and De Michele')

    % Larva 2 - Development

figure('Name', 'Larva 2')
hold on

boxchart(Temp_L2, L2_Rates)
plot(TempStore_L2, BriStore_L2, '-.', 'LineWidth', 1.25)
plot(TempStore_L2, LogStore_L2, '--', 'LineWidth', 1.25)
plot(TempStore_L2, LacOneStore_L2, ':', 'LineWidth', 1.25)
plot(TempStore_L2, LacTwoStore_L2, '-', 'LineWidth', 1.25)
plot(TempStore_L2, SDMStore_L2, 'LineWidth', 1.25)

title('L2 to L3 development')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 40]);
ylim([0 1.6]);
legend('Experimental data', 'Briere', 'Logan', ' Lactin-1', 'Lactin-2', ...
       'Sharpe and De Michele')

    % Larva 3 - Development

figure('Name', 'Larva 3')
hold on

boxchart(Temp_L3, L3_Rates)
plot(TempStore_L3, BriStore_L3, '-.', 'LineWidth', 1.25)
plot(TempStore_L3, LogStore_L3, '--', 'LineWidth', 1.25)
plot(TempStore_L3, LacOneStore_L3, ':', 'LineWidth', 1.25)
plot(TempStore_L3, LacTwoStore_L3, '-', 'LineWidth', 1.25)
plot(TempStore_L3, SDMStore_L3, 'LineWidth', 1.25)

title('L3 to pupa development')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 40]);
ylim([0 1]);
legend('Experimental data', 'Briere', 'Logan', ' Lactin-1', 'Lactin-2', ...
       'Sharpe and De Michele')

    % Pupa - Development

figure('Name', 'Pupa')
hold on

boxchart(Temp_P, P_Rates)
plot(TempStore_P, BriStore_P, '-.', 'LineWidth', 1.25)
plot(TempStore_P, LogStore_P, '--', 'LineWidth', 1.25)
plot(TempStore_P, LacOneStore_P, ':', 'LineWidth', 1.25)
plot(TempStore_P, LacTwoStore_P, '-', 'LineWidth', 1.25)
plot(TempStore_P, SDMStore_P, 'LineWidth', 1.25)

title('Pupa to adult development')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 40]);
ylim([0 0.3]);
legend('Experimental data', 'Briere', 'Logan', ' Lactin-1', 'Lactin-2', ...
       'Sharpe and De Michele')

   % Adult males survival

figure('Name', 'Adult males survival')
hold on
 
boxchart(Temp_AM, AM_Rates)
plot(TempStore_AM, WangStore_AM, '-.', 'LineWidth', 1.25)
plot(TempStore_AM, BriStore_AM, '-.', 'LineWidth', 1.25)
plot(TempStore_AM, LogStore_AM, '--', 'LineWidth', 1.25)
plot(TempStore_AM, LacOneStore_AM, ':', 'LineWidth', 1.25)
plot(TempStore_AM, LacTwoStore_AM, '-', 'LineWidth', 1.25)
plot(TempStore_AM, SDMStore_AM, 'LineWidth', 1.25)

title('Adult males survival')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 39]);
ylim([0 0.18]);
legend('Experimental data', 'Wang', 'Briere', 'Logan', ' Lactin-1', 'Lactin-2', ...
       'Sharpe and De Michele')

    % Adult females survival
 
figure('Name', 'Adult females survival')
hold on
 
boxchart(Temp_AF, AF_Rates)
plot(TempStore_AF, WangStore_AF, '-.', 'LineWidth', 1.25)
plot(TempStore_AF, BriStore_AF, '-.', 'LineWidth', 1.25)
plot(TempStore_AF, LogStore_AF, '--', 'LineWidth', 1.25)
plot(TempStore_AF, LacOneStore_AF, ':', 'LineWidth', 1.25)
plot(TempStore_AF, LacTwoStore_AF, '-', 'LineWidth', 1.25)
plot(TempStore_AF, SDMStore_AF, 'LineWidth', 1.25)

title('Adult females survival')
xlabel('Temperature (°C)')
ylabel('Development rate (1/days)')
xlim([0 39]);
ylim([0 0.18]);
legend('Experimental data', 'Wang', 'Briere', 'Logan', ' Lactin-1', ...
       'Lactin-2', 'Sharpe and De Michele')


    % Egg - Mortality

figure('Name', 'Egg')
hold on

scatter(Temp_Mortality, Mort_E, 'filled')
plot(TempStore_Wang, WangStore_E, '-.', 'LineWidth', 1.25)
plot(TempStore_KL, KLStore_E, 'LineWidth', 1.25)
plot(TempStore_ExpPol, ExpPolStore_E, 'LineWidth', 1.25)

title('Egg mortality')
xlabel('Temperature (°C)')
ylabel('Mortality rate (%)')
xlim([0 40]);
ylim([0 1.6]);
legend('Experimental data', 'Wang', 'Kim and Lee')

    % Larva 1 - Mortality

figure('Name', 'Larva 1')
hold on

scatter(Temp_Mortality, Mort_L1, 'filled')
plot(TempStore_Wang, WangStore_L1, '-.', 'LineWidth', 1.25)
plot(TempStore_KL, KLStore_L1, 'LineWidth', 1.25)
plot(TempStore_ExpPol, ExpPolStore_L1, 'LineWidth', 1.25)

title('L1 mortality')
xlabel('Temperature (°C)')
ylabel('Mortality rate (%)')
xlim([0 40]);
ylim([0 1.6]);
legend('Experimental data', 'Wang', 'Kim and Lee', 'Exponential')

    % Larva 2 - Mortality

figure('Name', 'Larva 2')
hold on

scatter(Temp_Mortality, Mort_L2, 'filled')
plot(TempStore_Wang, WangStore_L2, '-.', 'LineWidth', 1.25)
plot(TempStore_KL, KLStore_L2, 'LineWidth', 1.25)
plot(TempStore_ExpPol, ExpPolStore_L2, 'LineWidth', 1.25)

title('L2 mortality')
xlabel('Temperature (°C)')
ylabel('Mortality rate (%)')
xlim([0 40]);
ylim([0 1.6]);
legend('Experimental data', 'Wang', 'Kim and Lee', 'Exponential')

    % Larva 3 - Mortality

figure('Name', 'Larva 3')
hold on

scatter(Temp_Mortality, Mort_L3, 'filled')
plot(TempStore_Wang, WangStore_L3, '-.', 'LineWidth', 1.25)
plot(TempStore_KL, KLStore_L3, 'LineWidth', 1.25)
plot(TempStore_ExpPol, ExpPolStore_L3, 'LineWidth', 1.25)

title('L3 mortality')
xlabel('Temperature (°C)')
ylabel('Mortality rate (%)')
xlim([0 40]);
ylim([0 1.01]);
legend('Experimental data', 'Wang', 'Kim and Lee', 'Exponential')

    % Fertility rate

figure('Name', 'Total eggs per female')
hold on

scatter(Temp_Fertility, Tot_Egg_PerFemale, 'filled', 'Color', 'blue')
errorbar(Temp_Fertility, Tot_Egg_PerFemale, Err_Tot_Egg_PerFemale, ...
         'Color', 'blue')
plot(TempStore_Fert, FertStore, 'LineWidth', 1.25)

title('Total egg per female over temperature')
xlabel('Temperature (°C)')
ylabel('Totals eggs per female')
xlim([0 40]);
ylim([0 250]);
legend('Experimental data', 'Gaussian-like')

    % 3D scatterplot of egg production

figure('Name', 'Egg production over time and temperature')
hold on

scatter3(Temp_Fert_3D, TimeFert_3D, Egg_3D, 36, Temp_Fert_3D, 'filled')
colorbar
