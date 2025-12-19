
% Script where you can define initial values for the fits carried out in
% RunMe.m

% Created by Luca Rossini on 1 September 2024
% Last update: 3 September 2024

% E-mail: luca.rossini@ulb.be


%% Insert the initial conditions and the boundaries for the parameters of
%  each fit function

% Briere development rate function:

% Egg

InitCond_Bri_E = [1e-05 5.0 35.0 1.5];
LowBound_Bri_E = [1e-07 4.0 33.0 1.0];
UpBound_Bri_E = [0.01 10.0 40.0 4.0];

% L1

InitCond_Bri_L1 = [0.00012, 1.5, 30.5, 0.678];
LowBound_Bri_L1 = [0.00001 0.001 25 1.5];
UpBound_Bri_L1 = [0.0012 15.0 37.0 6.0];

% L2

InitCond_Bri_L2 = [0.00012, 1.5, 30.5, 0.678];
LowBound_Bri_L2 = [-0.0012 -15 -305 -6.78];
UpBound_Bri_L2 = [0.0012 15 305 6.78];

% L3

InitCond_Bri_L3 = [0.00012, 1.5, 30.5, 0.678];
LowBound_Bri_L3 = [-0.0012 -15 -305 -6.78];
UpBound_Bri_L3 = [0.0012 15 305 6.78];

% Pupa

InitCond_Bri_P = [1e-05 5.0 35.0 1.5];
LowBound_Bri_P = [1e-07 4.0 33.0 1.0];
UpBound_Bri_P = [0.01 10.0 40.0 4.0];

% Adult males

InitCond_Bri_AM = [1e-05 5.0 35.0 1.5];
LowBound_Bri_AM = [1e-07 0.01 33.0 0.5];
UpBound_Bri_AM = [0.01 10.0 40.0 4.0];

% Adult females

InitCond_Bri_AF = [1e-05 5.0 35.0 1.5];
LowBound_Bri_AF = [1e-07 0.01 33.0 0.5];
UpBound_Bri_AF = [0.01 10.0 40.0 4.0];

% Logan development rate function:

% Egg
        
InitCond_Log_E = [0.0679 0.167 32.91 4.463];
LowBound_Log_E = [-0.68 -1.67 -329.1 -44.63];
UpBound_Log_E = [0.68 1.67 329.1 44.63];

% L1

InitCond_Log_L1 = [0.218 0.187 32.667 3.779];
LowBound_Log_L1 = [-2.18 -1.87 -326.67 -37.79];
UpBound_Log_L1 = [2.18 1.87 326.67 37.79];

% L2
        
InitCond_Log_L2 = [0.218 0.187 32.667 3.779];
LowBound_Log_L2 = [-2.18 -1.87 -326.67 -37.79];
UpBound_Log_L2 = [2.18 1.87 326.67 37.79];

% L3
        
InitCond_Log_L3 = [0.218 0.187 32.667 3.779];
LowBound_Log_L3 = [-2.18 -1.87 -326.67 -37.79];
UpBound_Log_L3 = [2.18 1.87 326.67 37.79];

% Pupa
        
InitCond_Log_P = [0.0679 0.167 32.91 4.463];
LowBound_Log_P = [-0.68 -1.67 -329.1 -44.63];
UpBound_Log_P = [0.68 1.67 329.1 44.63];

% Adult males
        
InitCond_Log_AM = [0.001 0.067 32.91 4.463];
LowBound_Log_AM = [0.00001 0.0007 26.1 1.5];
UpBound_Log_AM = [0.01 0.9 37.1 6];

% Adult females
        
InitCond_Log_AF = [0.0679 0.067 32.91 4.463];
LowBound_Log_AF = [0.00001 0.0007 26.1 1.5];
UpBound_Log_AF = [0.01 0.9 37.1 6];


% Lactin-1 development rate function:

% Egg
        
InitCond_LacOne_E = [0.1422 33.976 6.8691];
LowBound_LacOne_E = [-1.40 -330 -69];
UpBound_LacOne_E = [1.40 330 69];

% L1
        
InitCond_LacOne_L1 = [0.1483 35.1585 6.6560];
LowBound_LacOne_L1 = [-1.40 -330 -69];
UpBound_LacOne_L1 = [1.40 330 69];

% L2
        
InitCond_LacOne_L2 = [0.1449 34.8657 6.8120];
LowBound_LacOne_L2 = [-1.40 -330 -69];
UpBound_LacOne_L2 = [1.40 330 69];

% L3
        
InitCond_LacOne_L3 = [0.1435 32.9023 6.8695];
LowBound_LacOne_L3 = [-1.40 -330 -69];
UpBound_LacOne_L3 = [1.40 330 69];

% Pupa
        
InitCond_LacOne_P = [0.1429 32.6309 6.9643];
LowBound_LacOne_P = [-1.40 -330 -69];
UpBound_LacOne_P = [1.40 330 69];

% Adult males
        
InitCond_LacOne_AM = [0.14 32.63 6.9];
LowBound_LacOne_AM = [0.0001 25 2];
UpBound_LacOne_AM = [1.40 35 12];

% Adult females
        
InitCond_LacOne_AF = [0.14 32.63 6.9];
LowBound_LacOne_AF = [0.0001 25 2];
UpBound_LacOne_AF = [1.40 35 12];


% Lactin-2 development rate function:

% Egg
        
InitCond_LacTwo_E = [0.1422 33.976 6.8691 -1.1001];
LowBound_LacTwo_E = [-0.361 -350.27 -31.90 -12.29];
UpBound_LacTwo_E = [0.361 350.27 31.90 12.29];

% L1
        
InitCond_LacTwo_L1 = [0.0260 33.8062 1.3973 -1.1001];
LowBound_LacTwo_L1 = [0.001 25.0 0.001 -12.29];
UpBound_LacTwo_L1 = [0.361 45.27 31.90 12.29];

% L2
        
InitCond_LacTwo_L2 = [0.0260 34.8169 2.0344 0.0001];
LowBound_LacTwo_L2 = [0.0001 25.0 0.001 0.0];
UpBound_LacTwo_L2 = [0.361 40.27 31.90 12.29];

% L3
        
InitCond_LacTwo_L3 = [0.0222 33.2772 1.46726 -1.14923];
LowBound_LacTwo_L3 = [-0.361 -350.27 -31.90 -12.29];
UpBound_LacTwo_L3 = [0.361 350.27 31.90 12.29];

% Pupa
        
InitCond_LacTwo_P = [0.0122 34.165614 2.1929 -1.10279];
LowBound_LacTwo_P = [-0.361 -350.27 -31.90 -12.29];
UpBound_LacTwo_P = [0.361 350.27 31.90 12.29];

% Adult males
        
InitCond_LacTwo_AM = [0.0001 32.63 6.9 -2];
LowBound_LacTwo_AM = [0.00001 25.0 2 -15];
UpBound_LacTwo_AM = [0.9 35.0 18 5];

% Adult females
        
InitCond_LacTwo_AF = [0.0001 32.63 6.9 -2];
LowBound_LacTwo_AF = [0.00001 25.0 2 -15];
UpBound_LacTwo_AF = [0.9 35.0 18 5];


% Sharpe and De Michele development rate function:

% Egg
        
InitCond_SDM_E = [4.1 -115.52 20.01 282.45 6.37 -134.08];
LowBound_SDM_E = [-40.71 -1151.50 -199.90 -2823.20 -63.70 -1338.10];
UpBound_SDM_E = [40.71 1151.5 199.90 2823.20 63.70 1338.10];

% L1
        
InitCond_SDM_L1 = [3.92 -118.10 20.19 299.89 6.64 -134.143];
LowBound_SDM_L1 = [-40.71 -1151.50 -199.90 -2823.2 -63.70 -1338.10];
UpBound_SDM_L1 = [40.71 1151.5 199.90 2823.20 63.70 1338.10];

% L2
        
InitCond_SDM_L2 = [4.32 -122.03 21.24 307.50 6.86 -146.32];
LowBound_SDM_L2 = [-40.71 -1151.50 -199.90 -2823.2 -63.70 -1338.10];
UpBound_SDM_L2 = [40.71 1151.5 199.90 2823.20 63.70 1338.10];

% L3
        
InitCond_SDM_L3 = [4.32 -122.03 21.24 307.60 6.86 -146.32];
LowBound_SDM_L3 = [-40.71 -1151.50 -199.90 -2823.2 -63.70 -1338.10];
UpBound_SDM_L3 = [40.71 1151.5 199.90 2823.20 63.70 1338.10];

% Pupa
        
InitCond_SDM_P = [7.15 -101.10 52.27 1091.39 11.49 -110.92];
LowBound_SDM_P = [-40.71 -1151.50 -199.90 -2823.2 -63.70 -1338.10];
UpBound_SDM_P = [40.71 1151.5 199.90 2823.20 63.70 1338.10];

% Adult males
        
InitCond_SDM_AM = [7.0 -70 14 8 10 -105];
LowBound_SDM_AM = [2 -130 9 4 5 -120];
UpBound_SDM_AM = [15 -50 25 12 12 -70];

% Adult females
        
InitCond_SDM_AF = [7.0 -70 14 8 10 -105];
LowBound_SDM_AF = [2 -130 9 4 5 -120];
UpBound_SDM_AF = [15 -50 25 12 12 -70];


% Wang mortality rate function

% Egg

InitCond_Wang_E = [1 1 1 1 3];
LowBound_Wang_E = [-3 -3 -3 -3 2];
UpBound_Wang_E = [100 100 100 100 10];

% L1

InitCond_Wang_L1 = [1 1 1 1 1];
LowBound_Wang_L1 = [-3 -3 -3 -3 -3];
UpBound_Wang_L1 = [100 100 100 100 100];

% L2

InitCond_Wang_L2 = [1 1 1 1 1];
LowBound_Wang_L2 = [-3 -3 -3 -3 -3];
UpBound_Wang_L2 = [100 100 100 100 100];

% L3

InitCond_Wang_L3 = [0.00003 -0.001 0.1 -1 5];
LowBound_Wang_L3 = [0 -1 0 -2 3];
UpBound_Wang_L3 = [1 2 3 0 7];

% Adult males

InitCond_Wang_AM = [1 1 1 1 1];
LowBound_Wang_AM = [-3 -3 -3 -3 -3];
UpBound_Wang_AM = [100 100 100 100 100];

% Adult females

InitCond_Wang_AF = [1 1 1 1 1];
LowBound_Wang_AF = [-3 -3 -3 -3 -3];
UpBound_Wang_AF = [100 100 100 100 100];


% Kim and Lee mortality rate function

% Egg

InitCond_KL_E = [0.9 23 8];
LowBound_KL_E = [0.7 17 5];
UpBound_KL_E = [1 25 11];

% L1

InitCond_KL_L1 = [1 1 1];
LowBound_KL_L1 = [-3 -3 -3 -3];
UpBound_KL_L1 = [100 100 100 100];

% L2

InitCond_KL_L2 = [0.8 16 9];
LowBound_KL_L2 = [0.7 14 6];
UpBound_KL_L2 = [1 22 12];

% L3

InitCond_KL_L3 = [0.8 16 9];
LowBound_KL_L3 = [0.7 14 6];
UpBound_KL_L3 = [1 23 12];


% Exponential poly mortality rate function

% Egg

InitCond_ExpPol_E = [0.01 0.3 1];
LowBound_ExpPol_E = [0 -Inf -Inf];
UpBound_ExpPol_E = [0.1 Inf Inf];

% L1

InitCond_ExpPol_L1 = [0.01 0.3 1];
LowBound_ExpPol_L1 = [0 -Inf -Inf];
UpBound_ExpPol_L1 = [0.1 Inf Inf];

% L2

InitCond_ExpPol_L2 = [0.01 0.3 1];
LowBound_ExpPol_L2 = [0 -Inf -Inf];
UpBound_ExpPol_L2 = [0.1 Inf Inf];

% L3

InitCond_ExpPol_L3 = [0.01 0.3 1];
LowBound_ExpPol_L3 = [0 -Inf -Inf];
UpBound_ExpPol_L3 = [0.1 Inf Inf];


% Egg production per female over temperature

InitCond_Fert = [33000., 80., 50., 6., 23.];
LowBound_Fert = [900., 10., 20., 4., 19.];
UpBound_Fert = [37000., 100., 70., 15., 35.];