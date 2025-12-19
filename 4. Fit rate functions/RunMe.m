
% Script for LifeTables analysis - Drosophila suzukii
% Created by Luca Rossini on 11 August 2024
% Last update 3 September 2024
% e-mail luca.rossini@ulb.be

clc
close all
clear

% Start to calculate the simulation time

tic

%% Acquisition of the data from 'DataAcquisition.m'

run('DataAcquisition.m')


%% Acquisition of the initial and boundary conditions for the fitting from 
% 'InitialParameters.m'

run('InitialParameters.m')


%% Fit EGG - Development

Stage_E = 'Egg';

% Briere : it calls the specific Briere fitting function

BriRes_E = Functions.BriFit(Temp_E, E_Rates, InitCond_Bri_E, LowBound_Bri_E, ...
                  UpBound_Bri_E, Stage_E);

% Logan : it calls the specific Logan fitting function

LogRes_E = Functions.LogFit(Temp_E, E_Rates, InitCond_Log_E, LowBound_Log_E, ...
                  UpBound_Log_E, Stage_E);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_E = Functions.LacOneFit(Temp_E, E_Rates, InitCond_LacOne_E, ...
                  LowBound_LacOne_E, UpBound_LacOne_E, Stage_E);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_E = Functions.LacTwoFit(Temp_E, E_Rates, InitCond_LacTwo_E, ...
                  LowBound_LacTwo_E, UpBound_LacTwo_E, Stage_E);

% Shape and De Michele : it calls the specific SDM fitting function

SDMRes_E = Functions.SDMFit(Temp_E, E_Rates, InitCond_SDM_E, LowBound_SDM_E, ...
                  UpBound_SDM_E, Stage_E);


%% Fit L1 - Development

Stage_L1 = 'Larva 1';

% Briere : it calls the specific Briere fitting function

BriRes_L1 = Functions.BriFit(Temp_L1, L1_Rates, InitCond_Bri_L1, ...
                             LowBound_Bri_L1, UpBound_Bri_L1, Stage_L1);

% Logan : it calls the specific Logan fitting function

LogRes_L1 = Functions.LogFit(Temp_L1, L1_Rates, InitCond_Log_L1, ...
                             LowBound_Log_L1, UpBound_Log_L1, Stage_L1);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_L1 = Functions.LacOneFit(Temp_L1, L1_Rates, InitCond_LacOne_L1, ...
                  LowBound_LacOne_L1, UpBound_LacOne_L1, Stage_L1);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_L1 = Functions.LacTwoFit(Temp_L1, L1_Rates, InitCond_LacTwo_L1, ...
                  LowBound_LacTwo_L1, UpBound_LacTwo_L1, Stage_L1);

% Shape and De Michele : it calls the specific SDM fitting function

SDMRes_L1 = Functions.SDMFit(Temp_L1, L1_Rates, InitCond_SDM_L1, ...
                             LowBound_SDM_L1, UpBound_SDM_L1, Stage_L1);


%% Fit L2 - Development

Stage_L2 = 'Larva 2';

% Briere : it calls the specific Briere fitting function

BriRes_L2 = Functions.BriFit(Temp_L2, L2_Rates, InitCond_Bri_L2, ...
                             LowBound_Bri_L2, UpBound_Bri_L2, Stage_L2);

% Logan : it calls the specific Logan fitting function

LogRes_L2 = Functions.LogFit(Temp_L2, L2_Rates, InitCond_Log_L2, ...
                             LowBound_Log_L2, UpBound_Log_L2, Stage_L2);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_L2 = Functions.LacOneFit(Temp_L2, L2_Rates, InitCond_LacOne_L2, ...
                  LowBound_LacOne_L2, UpBound_LacOne_L2, Stage_L2);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_L2 = Functions.LacTwoFit(Temp_L2, L2_Rates, InitCond_LacTwo_L2, ...
                  LowBound_LacTwo_L2, UpBound_LacTwo_L2, Stage_L2);

% Shape and De Michele : it calls the specific SDM fitting function

SDMRes_L2 = Functions.SDMFit(Temp_L2, L2_Rates, InitCond_SDM_L2, ...
                             LowBound_SDM_L2, UpBound_SDM_L2, Stage_L2);


%% Fit L3 - Development

Stage_L3 = 'Larva 3';

% Briere : it calls the specific Briere fitting function

BriRes_L3 = Functions.BriFit(Temp_L3, L3_Rates, InitCond_Bri_L3, ...
                             LowBound_Bri_L3, UpBound_Bri_L3, Stage_L3);

% Logan : it calls the specific Logan fitting function

LogRes_L3 = Functions.LogFit(Temp_L3, L3_Rates, InitCond_Log_L3, ...
                             LowBound_Log_L3, UpBound_Log_L3, Stage_L3);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_L3 = Functions.LacOneFit(Temp_L3, L3_Rates, InitCond_LacOne_L3, ...
                  LowBound_LacOne_L3, UpBound_LacOne_L3, Stage_L3);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_L3 = Functions.LacTwoFit(Temp_L3, L3_Rates, InitCond_LacTwo_L3, ...
                  LowBound_LacTwo_L3, UpBound_LacTwo_L3, Stage_L3);

% Shape and De Michele : it calls the specific SDM fitting function

SDMRes_L3 = Functions.SDMFit(Temp_L3, L3_Rates, InitCond_SDM_L3, ...
                             LowBound_SDM_L3, UpBound_SDM_L3, Stage_L3);


%% Fit P - Development

Stage_P = 'Pupa';

% Briere : it calls the specific Briere fitting function

BriRes_P = Functions.BriFit(Temp_P, P_Rates, InitCond_Bri_P, ...
                            LowBound_Bri_P, UpBound_Bri_P, Stage_P);

% Logan : it calls the specific Logan fitting function

LogRes_P = Functions.LogFit(Temp_P, P_Rates, InitCond_Log_P, ...
                            LowBound_Log_P, UpBound_Log_P, Stage_P);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_P = Functions.LacOneFit(Temp_P, P_Rates, InitCond_LacOne_P, ...
                  LowBound_LacOne_P, UpBound_LacOne_P, Stage_P);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_P = Functions.LacTwoFit(Temp_P, P_Rates, InitCond_LacTwo_P, ...
                  LowBound_LacTwo_P, UpBound_LacTwo_P, Stage_P);

% Shape and De Michele : it calls the specific SDM fitting function

SDMRes_P = Functions.SDMFit(Temp_P, P_Rates, InitCond_SDM_P, ...
                             LowBound_SDM_P, UpBound_SDM_P, Stage_P);



%% Fit AM - Development
 
Stage_AM = 'Adult males';
 
% Wang: it calls the specific Wang fitting function
 
WangRes_AM = Functions.WangFit(Temp_AM, AM_Rates, InitCond_Wang_AM, ...
                              LowBound_Wang_AM, UpBound_Wang_AM, Stage_AM);

% Briere : it calls the specific Briere fitting function

BriRes_AM = Functions.BriFit(Temp_AM, AM_Rates, InitCond_Bri_AM, ...
                             LowBound_Bri_AM, UpBound_Bri_AM, Stage_AM);
 
% Logan : it calls the specific Logan fitting function
 
LogRes_AM = Functions.LogFit(Temp_AM, AM_Rates, InitCond_Log_AM, ...
                             LowBound_Log_AM, UpBound_Log_AM, Stage_AM);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_AM = Functions.LacOneFit(Temp_AM, AM_Rates, InitCond_LacOne_AM, ...
                  LowBound_LacOne_AM, UpBound_LacOne_AM, Stage_AM);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_AM = Functions.LacTwoFit(Temp_AM, AM_Rates, InitCond_LacTwo_AM, ...
                  LowBound_LacTwo_AM, UpBound_LacTwo_AM, Stage_AM);

% Sharpe and De Michele : it calls the specific SDM fitting function

SDMRes_AM = Functions.SDMFit(Temp_AM, AM_Rates, InitCond_SDM_AM, ...
                             LowBound_SDM_AM, UpBound_SDM_AM, Stage_AM);

 
%% Fit AF

Stage_AF = 'Adult females';

% Wang: it calls the specific Wang fitting function
 
WangRes_AF = Functions.WangFit(Temp_AF, AF_Rates, InitCond_Wang_AF, ...
                              LowBound_Wang_AF, UpBound_Wang_AF, Stage_AF);

% Briere : it calls the specific Briere fitting function

BriRes_AF = Functions.BriFit(Temp_AF, AF_Rates, InitCond_Bri_AF, ...
                            LowBound_Bri_AF, UpBound_Bri_AF, Stage_AF);

% Logan : it calls the specific Logan fitting function

LogRes_AF = Functions.LogFit(Temp_AF, AF_Rates, InitCond_Log_AF, ...
                             LowBound_Log_AF, UpBound_Log_AF, Stage_AF);

% Lactin One : it calls the specific Lactin-1 fitting function

LacOneRes_AF = Functions.LacOneFit(Temp_AF, AF_Rates, InitCond_LacOne_AF, ...
                  LowBound_LacOne_AF, UpBound_LacOne_AF, Stage_AF);

% Lactin Two : it calls the specific Lactin-2 fitting function

LacTwoRes_AF = Functions.LacTwoFit(Temp_AF, AF_Rates, InitCond_LacTwo_AF, ...
                  LowBound_LacTwo_AF, UpBound_LacTwo_AF, Stage_AF);

% Sharpe and De Michele : it calls the specific SDM fitting function

SDMRes_AF = Functions.SDMFit(Temp_AF, AF_Rates, InitCond_SDM_AF, ...
                             LowBound_SDM_AF, UpBound_SDM_AF, Stage_AF);


%% Fitting the mortality rates

% Egg

% Wang

WangRes_E = Functions.WangFit(Temp_Mortality, Mort_E, InitCond_Wang_E, ...
                            LowBound_Wang_E, UpBound_Wang_E, Stage_E);

% Kim and Lee

KLRes_E = Functions.KLFit(Temp_Mortality, Mort_E, InitCond_KL_E, ...
                            LowBound_KL_E, UpBound_KL_E, Stage_E);

% Exponential poly

ExpPolRes_E = Functions.Exp1Fit(Temp_Mortality, Mort_E, InitCond_ExpPol_E, ...
                                LowBound_ExpPol_E, UpBound_ExpPol_E, Stage_E);

% L1

% Wang

WangRes_L1 = Functions.WangFit(Temp_Mortality, Mort_L1, InitCond_Wang_L1, ...
                             LowBound_Wang_L1, UpBound_Wang_L1, Stage_L1);

% Kim and Lee

KLRes_L1 = Functions.KLFit(Temp_Mortality, Mort_L1, InitCond_KL_L1, ...
                           LowBound_KL_L1, UpBound_KL_L1, Stage_L1);

% Exponential poly

ExpPolRes_L1 = Functions.Exp1Fit(Temp_Mortality, Mort_L1, InitCond_ExpPol_L1, ...
                                LowBound_ExpPol_L1, UpBound_ExpPol_L1, Stage_L1);

% L2

% Wang

WangRes_L2 = Functions.WangFit(Temp_Mortality, Mort_L2, InitCond_Wang_L2, ...
                             LowBound_Wang_L2, UpBound_Wang_L2, Stage_L2);

% Kim and Lee

KLRes_L2 = Functions.KLFit(Temp_Mortality, Mort_L2, InitCond_KL_L2, ...
                           LowBound_KL_L2, UpBound_KL_L2, Stage_L2);

% Exponential poly

ExpPolRes_L2 = Functions.Exp1Fit(Temp_Mortality, Mort_L2, InitCond_ExpPol_L2, ...
                                LowBound_ExpPol_L2, UpBound_ExpPol_L2, Stage_L2);

% L3

% Wang

WangRes_L3 = Functions.WangFit(Temp_Mortality, Mort_L3, InitCond_Wang_L3, ...
                             LowBound_Wang_L3, UpBound_Wang_L3, Stage_L3);

% Kim and Lee

KLRes_L3 = Functions.KLFit(Temp_Mortality, Mort_L3, InitCond_KL_L3, ...
                           LowBound_KL_L3, UpBound_KL_L3, Stage_L3);

% Exponential poly

ExpPolRes_L3 = Functions.Exp1Fit(Temp_Mortality, Mort_L3, InitCond_ExpPol_L3, ...
                                LowBound_ExpPol_L3, UpBound_ExpPol_L3, Stage_L3);


%% Fitting the egg production

% Gaussian-like: it calls the specific Wang fitting function

FertRes = Functions.FertFit(Temp_Fertility, Tot_Egg_PerFemale, ...
                    InitCond_Fert, LowBound_Fert, UpBound_Fert, ...
                    'Egg production');


%% Making plots by calling the dedicated scripts

run("Plots.m")

% End the calculation of the simulation time

toc