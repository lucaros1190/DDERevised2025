
% Script that absorb the data from the file 'LifeTablesDataset.xlsx'
% This file support the main script that fit the most common rate functions

% Created by Luca Rossini on 9 August 2024
% Last update: 9 August 2024

% E-mail: luca.rossini@ulb.be


%% Section 1: Definition of the dataset file parameters

% As a first step we have to indicate the name of the file containing the 
% dataset and of the specific tabs where the data we need to absorb are 
% stored. 

    % Name of the Excel file: change in case you change the name of the
    % Excel file!!!!

DataFile = 'LifeTablesDataset.xlsx';

    % Name of the tab within the file

IndLifeHistory = 'Individual-LifeHistory';


    % Import the names of the stages and of the experimental conditions 
    % (the temperatures of rearing). They will be used to create names for 
    % all the variables. 

StageNames = readcell(DataFile, 'Sheet', 'DatasetOverview', 'Range', ...
                      'C40:C49'); 

LenStageNames = length(StageNames);

    % Remove the spaces from the names of the conditions

    % NOTE: select only the cells that correspond to the experimental 
    % conditions. If a cell is empty, you will get an error.

CondNames = readcell(DataFile, 'Sheet', 'DatasetOverview', 'Range', ...
                     'C13:C26');
CondNames = erase(CondNames, ' ');
StageNames = replace(StageNames, '-', 'NA');

    % Convert numeric values to strings (will only be used if the name
    % consists only of numbers)

for i = 1:numel(CondNames)
    if isnumeric(CondNames{i})
        CondNames{i} = num2str(CondNames{i});
    end
end

    % Extract the data from the 'PopulationDynamics' tab - Preimaginal
    % stages

dataTable_Juvenile = readmatrix(DataFile, 'Sheet', IndLifeHistory, ...
                                'Range', 'C22:GL1021');

    % Extract the data from the 'PopulationDynamics' tab - Adult survival

dataTable_Adults = readmatrix(DataFile, 'Sheet', IndLifeHistory, ...
                                'Range', 'C1027:GL2026');


%% Section 2: Creation of the variables names and data acquisition 

    % Only juvenile stages

RawDataArray = cell(LenStageNames * length(CondNames), 2);

for k = 1:numel(CondNames)
    for i = 1:numel(StageNames)

            % Form variable name

            varName_Time = [StageNames{i} '_' CondNames{k} '_Time'];
            varName_Rate = [StageNames{i} '_' CondNames{k} '_Rate'];

            % Assign data to the dynamically created variable

            if strncmpi(StageNames{i}, 'A', 1) == 1
                eval([varName_Time ' = dataTable_Adults(:, i+(k-1) * 13);']);
            else
                eval([varName_Time ' = dataTable_Juvenile(:, i+(k-1) * 13);']);
            end 
            
            % Remove zeros and NaN

            eval([varName_Time ' = ' varName_Time '(~isnan(' varName_Time '));']);
            eval([varName_Time ' = nonzeros(' varName_Time ');']);
            
            % Convert times in rates

            eval([varName_Rate ' = 1 ./ ' varName_Time ';']);
    end
end


%% Preparation of the dataset per stage, ready for the fitting

% Egg stage

    % Temperatures

T6C_E = ones(length(E_6C_Time), 1) * 6;
T9C_E = ones(length(E_9C_Time), 1) * 9;
T13C_E = ones(length(E_13C_Time), 1) * 13;
T18C_E = ones(length(E_18C_Time), 1) * 18;
T20C_E = ones(length(E_20C_Time), 1) * 20;
T24C_E = ones(length(E_24C_Time), 1) * 24;
T25C_E = ones(length(E_25C_Time), 1) * 25;
T26C_E = ones(length(E_26C_Time), 1) * 26;
T27C_E = ones(length(E_27C_Time), 1) * 27;
T28C_E = ones(length(E_28C_Time), 1) * 28;
T29C_E = ones(length(E_29C_Time), 1) * 29;
T31C_E = ones(length(E_31C_Time), 1) * 31;
T32C_E = ones(length(E_32C_Time), 1) * 32;
T33C_E = ones(length(E_33C_Time), 1) * 33;

Temp_E = [T6C_E; T9C_E; T13C_E; T18C_E; T20C_E; T24C_E; T25C_E; T26C_E; 
           T27C_E; T28C_E; T29C_E; T31C_E; T32C_E; T33C_E];

    % Egg array

E_Rates = [E_6C_Rate; E_9C_Rate; E_13C_Rate; E_18C_Rate; E_20C_Rate;
             E_24C_Rate; E_25C_Rate; E_26C_Rate; E_27C_Rate; E_28C_Rate;
             E_29C_Rate; E_31C_Rate; E_32C_Rate; E_33C_Rate];


% L1 stage

    % Temperatures

T6C_L1 = ones(length(L1_6C_Time), 1) * 6;
T9C_L1 = ones(length(L1_9C_Time), 1) * 9;
T13C_L1 = ones(length(L1_13C_Time), 1) * 13;
T18C_L1 = ones(length(L1_18C_Time), 1) * 18;
T20C_L1 = ones(length(L1_20C_Time), 1) * 20;
T24C_L1 = ones(length(L1_24C_Time), 1) * 24;
T25C_L1 = ones(length(L1_25C_Time), 1) * 25;
T26C_L1 = ones(length(L1_26C_Time), 1) * 26;
T27C_L1 = ones(length(L1_27C_Time), 1) * 27;
T28C_L1 = ones(length(L1_28C_Time), 1) * 28;
T29C_L1 = ones(length(L1_29C_Time), 1) * 29;
T31C_L1 = ones(length(L1_31C_Time), 1) * 31;
T32C_L1 = ones(length(L1_32C_Time), 1) * 32;
T33C_L1 = ones(length(L1_33C_Time), 1) * 33;

Temp_L1 = [T6C_L1; T9C_L1; T13C_L1; T18C_L1; T20C_L1; T24C_L1; T25C_L1; 
          T26C_L1; T27C_L1; T28C_L1; T29C_L1; T31C_L1; T32C_L1; T33C_L1];

    % L1 array

L1_Rates = [L1_6C_Rate; L1_9C_Rate; L1_13C_Rate; L1_18C_Rate; L1_20C_Rate;
            L1_24C_Rate; L1_25C_Rate; L1_26C_Rate; L1_27C_Rate; L1_28C_Rate;
            L1_29C_Rate; L1_31C_Rate; L1_32C_Rate; L1_33C_Rate];


% L2 stage

    % Temperatures

T6C_L2 = ones(length(L2_6C_Time), 1) * 6;
T9C_L2 = ones(length(L2_9C_Time), 1) * 9;
T13C_L2 = ones(length(L2_13C_Time), 1) * 13;
T18C_L2 = ones(length(L2_18C_Time), 1) * 18;
T20C_L2 = ones(length(L2_20C_Time), 1) * 20;
T24C_L2 = ones(length(L2_24C_Time), 1) * 24;
T25C_L2 = ones(length(L2_25C_Time), 1) * 25;
T26C_L2 = ones(length(L2_26C_Time), 1) * 26;
T27C_L2 = ones(length(L2_27C_Time), 1) * 27;
T28C_L2 = ones(length(L2_28C_Time), 1) * 28;
T29C_L2 = ones(length(L2_29C_Time), 1) * 29;
T31C_L2 = ones(length(L2_31C_Time), 1) * 31;
T32C_L2 = ones(length(L2_32C_Time), 1) * 32;
T33C_L2 = ones(length(L2_33C_Time), 1) * 33;

Temp_L2 = [T6C_L2; T9C_L2; T13C_L2; T18C_L2; T20C_L2; T24C_L2; T25C_L2; 
          T26C_L2; T27C_L2; T28C_L2; T29C_L2; T31C_L2; T32C_L2; T33C_L2];

    % L2 array

L2_Rates = [L2_6C_Rate; L2_9C_Rate; L2_13C_Rate; L2_18C_Rate; L2_20C_Rate;
            L2_24C_Rate; L2_25C_Rate; L2_26C_Rate; L2_27C_Rate; L2_28C_Rate;
            L2_29C_Rate; L2_31C_Rate; L2_32C_Rate; L2_33C_Rate];


% L3 stage

    % Temperatures

T6C_L3 = ones(length(L3_6C_Time), 1) * 6;
T9C_L3 = ones(length(L3_9C_Time), 1) * 9;
T13C_L3 = ones(length(L3_13C_Time), 1) * 13;
T18C_L3 = ones(length(L3_18C_Time), 1) * 18;
T20C_L3 = ones(length(L3_20C_Time), 1) * 20;
T24C_L3 = ones(length(L3_24C_Time), 1) * 24;
T25C_L3 = ones(length(L3_25C_Time), 1) * 25;
T26C_L3 = ones(length(L3_26C_Time), 1) * 26;
T27C_L3 = ones(length(L3_27C_Time), 1) * 27;
T28C_L3 = ones(length(L3_28C_Time), 1) * 28;
T29C_L3 = ones(length(L3_29C_Time), 1) * 29;
T31C_L3 = ones(length(L3_31C_Time), 1) * 31;
T32C_L3 = ones(length(L3_32C_Time), 1) * 32;
T33C_L3 = ones(length(L3_33C_Time), 1) * 33;

Temp_L3 = [T6C_L3; T9C_L3; T13C_L3; T18C_L3; T20C_L3; T24C_L3; T25C_L3; 
          T26C_L3; T27C_L3; T28C_L3; T29C_L3; T31C_L3; T32C_L3; T33C_L3];

    % L3 array

L3_Rates = [L3_6C_Rate; L3_9C_Rate; L3_13C_Rate; L3_18C_Rate; L3_20C_Rate;
            L3_24C_Rate; L3_25C_Rate; L3_26C_Rate; L3_27C_Rate; L3_28C_Rate;
            L3_29C_Rate; L3_31C_Rate; L3_32C_Rate; L3_33C_Rate];


% Pupa stage

    % Temperatures

T6C_P = ones(length(P_6C_Time), 1) * 6;
T9C_P = ones(length(P_9C_Time), 1) * 9;
T13C_P = ones(length(P_13C_Time), 1) * 13;
T18C_P = ones(length(P_18C_Time), 1) * 18;
T20C_P = ones(length(P_20C_Time), 1) * 20;
T24C_P = ones(length(P_24C_Time), 1) * 24;
T25C_P = ones(length(P_25C_Time), 1) * 25;
T26C_P = ones(length(P_26C_Time), 1) * 26;
T27C_P = ones(length(P_27C_Time), 1) * 27;
T28C_P = ones(length(P_28C_Time), 1) * 28;
T29C_P = ones(length(P_29C_Time), 1) * 29;
T31C_P = ones(length(P_31C_Time), 1) * 31;
T32C_P = ones(length(P_32C_Time), 1) * 32;
T33C_P = ones(length(P_33C_Time), 1) * 33;

Temp_P = [T6C_P; T9C_P; T13C_P; T18C_P; T20C_P; T24C_P; T25C_P; 
          T26C_P; T27C_P; T28C_P; T29C_P; T31C_P; T32C_P; T33C_P];

    % Pupa array

P_Rates = [P_6C_Rate; P_9C_Rate; P_13C_Rate; P_18C_Rate; P_20C_Rate;
           P_24C_Rate; P_25C_Rate; P_26C_Rate; P_27C_Rate; P_28C_Rate;
           P_29C_Rate; P_31C_Rate; P_32C_Rate; P_33C_Rate];


% Adult males stage survival

    % Temperatures

T6C_AM = ones(length(AM_6C_Time), 1) * 6;
T9C_AM = ones(length(AM_9C_Time), 1) * 9;
T13C_AM = ones(length(AM_13C_Time), 1) * 13;
T18C_AM = ones(length(AM_18C_Time), 1) * 18;
T20C_AM = ones(length(AM_20C_Time), 1) * 20;
T24C_AM = ones(length(AM_24C_Time), 1) * 24;
T25C_AM = ones(length(AM_25C_Time), 1) * 25;
T26C_AM = ones(length(AM_26C_Time), 1) * 26;
T27C_AM = ones(length(AM_27C_Time), 1) * 27;
T28C_AM = ones(length(AM_28C_Time), 1) * 28;
T29C_AM = ones(length(AM_29C_Time), 1) * 29;
T31C_AM = ones(length(AM_31C_Time), 1) * 31;
T32C_AM = ones(length(AM_32C_Time), 1) * 32;
T33C_AM = ones(length(AM_33C_Time), 1) * 33;

Temp_AM = [T6C_AM; T9C_AM; T13C_AM; T18C_AM; T20C_AM; T24C_AM; T25C_AM; 
          T26C_AM; T27C_AM; T28C_AM; T29C_AM; T31C_AM; T32C_AM; T33C_AM];

    % Adult males array

AM_Rates = [AM_6C_Rate; AM_9C_Rate; AM_13C_Rate; AM_18C_Rate; AM_20C_Rate;
            AM_24C_Rate; AM_25C_Rate; AM_26C_Rate; AM_27C_Rate; AM_28C_Rate;
            AM_29C_Rate; AM_31C_Rate; AM_32C_Rate; AM_33C_Rate];


% Adult females stage survival

    % Temperatures

T6C_AF = ones(length(AF_6C_Time), 1) * 6;
T9C_AF = ones(length(AF_9C_Time), 1) * 9;
T13C_AF = ones(length(AF_13C_Time), 1) * 13;
T18C_AF = ones(length(AF_18C_Time), 1) * 18;
T20C_AF = ones(length(AF_20C_Time), 1) * 20;
T24C_AF = ones(length(AF_24C_Time), 1) * 24;
T25C_AF = ones(length(AF_25C_Time), 1) * 25;
T26C_AF = ones(length(AF_26C_Time), 1) * 26;
T27C_AF = ones(length(AF_27C_Time), 1) * 27;
T28C_AF = ones(length(AF_28C_Time), 1) * 28;
T29C_AF = ones(length(AF_29C_Time), 1) * 29;
T31C_AF = ones(length(AF_31C_Time), 1) * 31;
T32C_AF = ones(length(AF_32C_Time), 1) * 32;
T33C_AF = ones(length(AF_33C_Time), 1) * 33;

Temp_AF = [T6C_AF; T9C_AF; T13C_AF; T18C_AF; T20C_AF; T24C_AF; T25C_AF; 
          T26C_AF; T27C_AF; T28C_AF; T29C_AF; T31C_AF; T32C_AF; T33C_AF];

    % Adult females array

AF_Rates = [AF_6C_Rate; AF_9C_Rate; AF_13C_Rate; AF_18C_Rate; AF_20C_Rate;
            AF_24C_Rate; AF_25C_Rate; AF_26C_Rate; AF_27C_Rate; AF_28C_Rate;
            AF_29C_Rate; AF_31C_Rate; AF_32C_Rate; AF_33C_Rate];


%% Absorbing mortality rates

DataFile_Mortality = 'Mortality.xlsx';

DataMort = readmatrix(DataFile_Mortality, 'Sheet', 'MortalityRates', 'Range', ...
                      'G2:K15');

Temp_Mortality = DataMort(:,1);

Mort_E = DataMort(:,2);
Mort_L1 = DataMort(:,3);
Mort_L2 = DataMort(:,4);
Mort_L3 = DataMort(:,5);


%% Absorbing egg production over temperature

DataFile_Fertility = 'Fertility.xlsx';

DataFert = readmatrix(DataFile_Fertility, 'Sheet', 'Analysis', 'Range', ...
                      'B4:D12');

Temp_Fertility = DataFert(:,1);
Tot_Egg_PerFemale = DataFert(:,2);
Err_Tot_Egg_PerFemale = DataFert(3);

% For 3D plot

DataFert_3D = readmatrix(DataFile_Fertility, 'Sheet', 'Analysis', 'Range', ...
                         'B35:D266');

Temp_Fert_3D = DataFert_3D(:,1);
TimeFert_3D = DataFert_3D(:,2);
Egg_3D = DataFert_3D(:,3);
