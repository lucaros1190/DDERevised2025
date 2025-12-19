
% Definition of the functions contained in the AnalisiStepMatlab.m file

% Created by Luca Rossini on 11 December 2023
% Last update 3 September 2024
% E-mail: luca.rossini@ulb.be


classdef Functions

    methods (Static) % Insert any function in between "methods" and "end".

    
    %% Briere fit function

    function BriPar = BriFit(Temperature, DevRates, InitCond, LowBound, ...
                             UpBound, Stage)


        % Definition of the mathematical expression
                % a = a
                % b = T_L
                % c = T_M
                % d = m

            BriereFun = ['abs(a * x * ( x - b ) * (( c - x )^( 1 / d )))'];

        % Initial conditions

            BriInitCond = InitCond;
            BriLowBound = LowBound;
            BriUpBound = UpBound;

        % Fit options

            BriFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                      'Lower', BriLowBound, 'Upper', ...
                                      BriUpBound, 'StartPoint', ...
                                      BriInitCond);

        % Fitting algorithm and solution storage

            [BriFit, Brigof] = fit(Temperature, DevRates, BriereFun, ...
                                   BriFitOpt);
    
        % Save the coeffiecient values for a, b, c, and d in a vector

            FitRes = coeffvalues(BriFit);

        % Get the confidence interval

            CI = confint(BriFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = BriFit(Temperature);
            Bins = 0 : (length(DevRates) - 1);
            
            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', DevRates, ...
                                 'Expected', FitData);
            
        % Print results on the console

            fprintf('\n Briere %s - Fit Results: \n\n', Stage)
            fprintf('a = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('T_L = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('T_M = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))
            fprintf('m = %f +/- %f\n', FitRes(4), ...
                                     abs((CI(7) - CI(8)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            Brigof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            BriPar = [FitRes(1) FitRes(2) FitRes(3) FitRes(4)];

        end


    %% Logan fit function

    function LogPar = LogFit(Temperature, DevRates, InitCond, LowBound, ...
                             UpBound, Stage)


        % Definition of the mathematical expression
                % a = psi
                % b = rho
                % c = T_M
                % d = Delta_T

            LoganFun = ['a * ( exp( b * x) - exp( b * c - ' ...
                        '( ( c - x) / d ) ) )'];

        % Initial conditions

            LogInitCond = InitCond;
            LogLowBound = LowBound;
            LogUpBound = UpBound;

        % Fit options

            LogFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                      'Lower', LogLowBound, 'Upper', ...
                                      LogUpBound, 'StartPoint', ...
                                      LogInitCond);

        % Fitting algorithm and solution storage

            [LogFit, Loggof] = fit(Temperature, DevRates, LoganFun, ...
                                   LogFitOpt);
    
        % Save the coeffiecient values for a, b, c, and d in a vector

            FitRes = coeffvalues(LogFit);

        % Get the confidence interval

            CI = confint(LogFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = LogFit(Temperature);
            Bins = 0 : (length(DevRates) - 1);

            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', DevRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf('\n Logan %s - Fit Results: \n\n', Stage)
            fprintf('psi = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('rho = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('T_M = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))
            fprintf('Delta_T = %f +/- %f\n', FitRes(4), ...
                                     abs((CI(7) - CI(8)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            Loggof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            LogPar = [FitRes(1) FitRes(2) FitRes(3) FitRes(4)];

        end


    %% Lactin-1 fit function

    function LacOnePar = LacOneFit(Temperature, DevRates, InitCond, ...
                                   LowBound, UpBound, Stage)


        % Definition of the mathematical expression
                % a = rho
                % b = T_M
                % c = Delta_T

            LacOneFun = ['exp( a * x) - exp( a * b - ( ( b - x ) / c ) )'];

        % Initial conditions

            LacOneInitCond = InitCond;
            LacOneLowBound = LowBound;
            LacOneUpBound = UpBound;

        % Fit options

            LacOneFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                      'Lower', LacOneLowBound, 'Upper', ...
                                      LacOneUpBound, 'StartPoint', ...
                                      LacOneInitCond);

        % Fitting algorithm and solution storage

            [LacOneFit, LacOnegof] = fit(Temperature, DevRates, LacOneFun, ...
                                         LacOneFitOpt);
    
        % Save the coeffiecient values for a, b, and c in a vector

            FitRes = coeffvalues(LacOneFit);

        % Get the confidence interval

            CI = confint(LacOneFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = LacOneFit(Temperature);
            Bins = 0 : (length(DevRates) - 1);

            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', DevRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf('\n Lactin-1 %s - Fit Results: \n\n', Stage)
            fprintf('rho = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('T_M = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('Delta_t = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            LacOnegof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            LacOnePar = [FitRes(1) FitRes(2) FitRes(3)];

        end


    %% Lactin-2 fit function

    function LacTwoPar = LacTwoFit(Temperature, DevRates, InitCond, ...
                                   LowBound, UpBound, Stage)


        % Definition of the mathematical expression
                % a = rho
                % b = T_M
                % c = Delta_T
                % d = lambda

            LacTwoFun = ['exp( a * x) - exp( a * b - ( ( b - x ) / c ) ) ' ...
                         ' + d'];

        % Initial conditions

            LacTwoInitCond = InitCond;
            LacTwoLowBound = LowBound;
            LacTwoUpBound = UpBound;

        % Fit options

            LacTwoFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                      'Lower', LacTwoLowBound, 'Upper', ...
                                      LacTwoUpBound, 'StartPoint', ...
                                      LacTwoInitCond);

        % Fitting algorithm and solution storage

            [LacTwoFit, LacTwogof] = fit(Temperature, DevRates, LacTwoFun, ...
                                         LacTwoFitOpt);
    
        % Save the coeffiecient values for a, b, c, and d in a vector

            FitRes = coeffvalues(LacTwoFit);

        % Get the confidence interval

            CI = confint(LacTwoFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = abs(LacTwoFit(Temperature));
            Bins = 0 : (length(DevRates) - 1);
            
            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', DevRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf('\n Lactin-2 %s - Fit Results: \n\n', Stage)
            fprintf('rho = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('T_M = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('Delta_t = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))
            fprintf('lambda = %f +/- %f\n', FitRes(4), ...
                                     abs((CI(7) - CI(8)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            LacTwogof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            LacTwoPar = [FitRes(1) FitRes(2) FitRes(3) FitRes(4)];

        end


    %% Wang (mortality) fit function

    function WangPar = WangFit(Temperature, MortRates, InitCond, LowBound, ...
                               UpBound, Stage)


        % Definition of the mathematical expression
                % a = a
                % b = b
                % c = c
                % d = d
                % e = e

            WangFun = ['a * x^4 + b * x^3 + c * x^2 + d * x + e'];

        % Initial conditions

            WangInitCond = InitCond;
            WangLowBound = LowBound;
            WangUpBound = UpBound;

        % Fit options

            WangFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                    'Lower', WangLowBound, 'Upper', ...
                                    WangUpBound, 'StartPoint', ...
                                    WangInitCond);

        % Fitting algorithm and solution storage

            [WangFit, Wanggof] = fit(Temperature, MortRates, WangFun, ...
                                         WangFitOpt);
    
        % Save the coeffiecient values for a, b, c, and d in a vector

            FitRes = coeffvalues(WangFit);

        % Get the confidence interval

            CI = confint(WangFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = abs(WangFit(Temperature));
            Bins = 0 : (length(MortRates) - 1);
            
            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', MortRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf('\n Wang %s - Fit Results: \n\n', Stage)
            fprintf('a = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('b = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('c = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))
            fprintf('d = %f +/- %f\n', FitRes(4), ...
                                     abs((CI(7) - CI(8)) / 3.92))
            fprintf('e = %f +/- %f\n', FitRes(5), ...
                                     abs((CI(9) - CI(10)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            Wanggof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            WangPar = [FitRes(1) FitRes(2) FitRes(3) FitRes(4) FitRes(5)];

        end


    %% Kim and Lee mortality fit function

    function KLPar = KLFit(Temperature, MortRates, InitCond, LowBound, ...
                           UpBound, Stage)


        % Definition of the mathematical expression
                % a = k
                % b = T_M
                % c = rho_T

            KLFun = ['1 - (a * exp(1 + ((b - x) / c) - exp((b - x) / c)))'];

        % Initial conditions

            KLInitCond = InitCond;
            KLLowBound = LowBound;
            KLUpBound = UpBound;

        % Fit options

            KLFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                  'Lower', KLLowBound, 'Upper', KLUpBound, ...
                                  'StartPoint', KLInitCond);

        % Fitting algorithm and solution storage

            [KLFit, KLgof] = fit(Temperature, MortRates, KLFun, KLFitOpt);
    
        % Save the coeffiecient values for a, b, and c in a vector

            FitRes = coeffvalues(KLFit);

        % Get the confidence interval

            CI = confint(KLFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = abs(KLFit(Temperature));
            Bins = 0 : (length(MortRates) - 1);
            
            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', MortRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf('\n Kim and Lee %s - Fit Results: \n\n', Stage)
            fprintf('k = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('T_M = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('rho_T = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            KLgof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            KLPar = [FitRes(1) FitRes(2) FitRes(3)];

        end


        %% Exponential-1 (mortality) fit function

        function Exp1Par = Exp1Fit(Temperature, MortRates, InitCond, LowBound, ...
                                   UpBound, Stage)


        % Definition of the mathematical expression
                % a = a
                % b = b
                % c = c

            Exp1Fun = ['exp(a * x^2 + b * x + c)'];

        % Initial conditions

            Exp1InitCond = InitCond;
            Exp1LowBound = LowBound;
            Exp1UpBound = UpBound;

        % Fit options

            Exp1FitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                     'Lower', Exp1LowBound, 'Upper', ...
                                     Exp1UpBound, 'StartPoint', ...
                                     Exp1InitCond);

        % Fitting algorithm and solution storage

            [Exp1Fit, Exp1gof] = fit(Temperature, MortRates, Exp1Fun, ...
                                     Exp1FitOpt);
    
        % Save the coeffiecient values for a, b, c, and d in a vector

            FitRes = coeffvalues(Exp1Fit);

        % Get the confidence interval

            CI = confint(Exp1Fit, 0.95);

        % Chi square test for the goodness of fit

            FitData = abs(Exp1Fit(Temperature));
            Bins = 0 : (length(MortRates) - 1);
            
            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', MortRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf(['\n Type Exponential-1 %s - Fit Results: \n' ...
                     'exp(a * x^2 + b * x + c) \n\n'], Stage)
            fprintf('a = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 3.92))
            fprintf('b = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 3.92))
            fprintf('c = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 3.92))

            fprintf('\n Goodness of fit:\n')
            Exp1gof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            Exp1Par = [FitRes(1) FitRes(2) FitRes(3)];

        end


    %% Gaussian-like (fertility) fit function

    function FertPar = FertFit(Temperature, FertRates, InitCond, LowBound, ...
                               UpBound, Stage)


       % Definition of the mathematical expression
               % a = alpha
               % b = gamma
               % c = lambda
               % d = tau
               % e = delta

           FertFun = ['a * ( ((b + 1) / (pi * (c ^ (2 * b + 2)))) * ' ...
                      '((c ^ 2) - ( ((x - e) ^ 2) + (d ^ 2)) ) ^ b )'];

       % Initial conditions

           FertInitCond = InitCond;
           FertLowBound = LowBound;
           FertUpBound = UpBound;

       % Fit options

           FertFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                   'Lower', FertLowBound, 'Upper', ...
                                   FertUpBound, 'StartPoint', FertInitCond);

       % Fitting algorithm and solution storage

           [FertFit, Fertgof] = fit(Temperature, FertRates, FertFun, ...
                                         FertFitOpt);
    
       % Save the coeffiecient values for a, b, c, and d in a vector

           FitRes = coeffvalues(FertFit);

       % Get the confidence interval

           CI = confint(FertFit, 0.95);

       % Chi square test for the goodness of fit

           FitData = abs(FertFit(Temperature));
           Bins = 0 : (length(FertRates) - 1);
            
           [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, 'Frequency', FertRates, ...
                                'Expected', FitData);

       % Print results on the console

           fprintf('\n Gaussian-like %s - Fit Results: \n\n', Stage)
           fprintf('alpha = %f +/- %f\n', FitRes(1), ...
                                    abs((CI(1) - CI(2)) / 3.92))
           fprintf('gamma = %f +/- %f\n', FitRes(2), ...
                                    abs((CI(3) - CI(4)) / 3.92))
           fprintf('lambda = %f +/- %f\n', FitRes(3), ...
                                    abs((CI(5) - CI(6)) / 3.92))
           fprintf('tau = %f +/- %f\n', FitRes(4), ...
                                    abs((CI(7) - CI(8)) / 3.92))
           fprintf('delta = %f +/- %f\n', FitRes(5), ...
                                    abs((CI(9) - CI(10)) / 3.92))

           fprintf('\n Goodness of fit:\n')
           Fertgof

           fprintf('\n Chi squared test:\n')
           h
           p
           st

           fprintf('\n')

           FertPar = [FitRes(1) FitRes(2) FitRes(3) FitRes(4) FitRes(5)];

       end


    %% Sharpe and De Michele fit function

    function SDMPar = SDMFit(Temperature, DevRates, InitCond, LowBound, ...
                             UpBound, Stage)


        % Definition of the mathematical expression

            SDMFun = ['(x * ( exp( A - (B / x ) ) ) ) / ' ...
                      '((1 + exp( C - ( D / x ) ) + exp( E - ( F / x ) ) ) )'];

        % Initial conditions

            SDMInitCond = InitCond;
            SDMLowBound = LowBound;
            SDMUpBound = UpBound;

        % Fit options

            SDMFitOpt = fitoptions('Method', 'NonlinearLeastSquares', ...
                                      'Lower', SDMLowBound, 'Upper', ...
                                      SDMUpBound, 'StartPoint', ...
                                      SDMInitCond);

        % Fitting algorithm and solution storage

            [SDMFit, SDMgof] = fit(Temperature, DevRates, SDMFun, SDMFitOpt);
    
        % Save the coeffiecient values in a vector

            FitRes = coeffvalues(SDMFit);

        % Get the confidence interval

            CI = confint(SDMFit, 0.95);

        % Chi square test for the goodness of fit

            FitData = SDMFit(Temperature);
            Bins = 0 : (length(DevRates) - 1);

            [h, p, st] = chi2gof(Bins, 'Ctrs', Bins, ...
                                 'Frequency', DevRates, ...
                                 'Expected', FitData);

        % Print results on the console

            fprintf('\n Sharpe and De Michele %s - Fit Results: \n\n', ...
                    Stage)
            fprintf('A = %f +/- %f\n', FitRes(1), ...
                                     abs((CI(1) - CI(2)) / 5.15))
            fprintf('B = %f +/- %f\n', FitRes(2), ...
                                     abs((CI(3) - CI(4)) / 5.15))
            fprintf('C = %f +/- %f\n', FitRes(3), ...
                                     abs((CI(5) - CI(6)) / 5.15))
            fprintf('D = %f +/- %f\n', FitRes(4), ...
                                     abs((CI(7) - CI(8)) / 5.15))
            fprintf('E = %f +/- %f\n', FitRes(5), ...
                                     abs((CI(9) - CI(10)) / 5.15))
            fprintf('F = %f +/- %f\n', FitRes(6), ...
                                     abs((CI(11) - CI(12)) / 5.15))

            fprintf('\n Goodness of fit:\n')
            SDMgof

            fprintf('\n Chi squared test:\n')
            h
            p
            st

            fprintf('\n')

            SDMPar = [FitRes(1) FitRes(2) FitRes(3) FitRes(4) FitRes(5) ...
                      FitRes(6)];

        end


    %% Evaluate the Brière for plotting purposes

    function Bri_TMP = Bri_Plot(T, BriPar)

        if T > BriPar(3)

            Bri_TMP = 0;
    
        else

            Bri_TMP = BriPar(1) * T * ( T - BriPar(2) ) * ...
                          (( BriPar(3) - T )^( 1 / BriPar(4) ));
        end

        end


    %% Evaluate the Logan for plotting purposes

    function Log_TMP = Log_Plot(T, LogPar)

        if T > LogPar(3)

            Log_TMP = 0;

        else

            Log_TMP = LogPar(1) * ( exp( LogPar(2) * T) - exp( LogPar(2) * ...
                      LogPar(3) - ( ( LogPar(3) - T) / LogPar(4) ) ) );

        end

    end


    %% Evaluate the Lactin-1 for plotting purposes

    function LacOne_TMP = LacOne_Plot(T, LacOnePar)

        if T > LacOnePar(2)

            LacOne_TMP = 0;

        else

            LacOne_TMP = exp( LacOnePar(1) * T) - ...
                        exp( LacOnePar(1) * LacOnePar(2) - ( ( LacOnePar(2) ...
                        - T ) / LacOnePar(3) ) );

    end

        end


    %% Evaluate the Lactin-2 for plotting purposes

    function LacTwo_TMP = LacTwo_Plot(T, LacTwoPar)
        
        if T > LacTwoPar(2)

            LacTwo_TMP = 0;

        else

            LacTwo_TMP = exp( LacTwoPar(1) * T) - exp( LacTwoPar(1) * ...
                         LacTwoPar(2) - ( ( LacTwoPar(2) - T ) / ...
                         LacTwoPar(3) ) ) + LacTwoPar(4);

        end

        if LacTwo_TMP < 0

            LacTwo_TMP = 0;

        else

            LacTwo_TMP = exp( LacTwoPar(1) * T) - exp( LacTwoPar(1) * ...
                         LacTwoPar(2) - ( ( LacTwoPar(2) - T ) / ...
                         LacTwoPar(3) ) ) + LacTwoPar(4);
    
        end

            end


    %% Evaluate the Sharpe and De Michele for plotting purposes

    function SDM_TMP = SDM_Plot(T, SDMPar)

        SDM_TMP = (T * ( exp( SDMPar(1) - (SDMPar(2) / T ) ) ) ) / ...
                 ((1 + exp( SDMPar(3) - ( SDMPar(4) / T ) ) + exp( SDMPar(5) - ...
                 ( SDMPar(6) / T ) ) ) );

        end


        %% Evaluate the Wang for plotting purposes

    function Wang_TMP = Wang_Plot(T, WangPar)

           Wang_TMP = WangPar(1) * T^4 + WangPar(2) * T^3 + WangPar(3) * ...
                       T^2 + WangPar(4) * T + WangPar(5);

        end

        
    %% Evaluate the Exponential 1 for plotting purposes

    function Exp1_TMP = Exp1_Plot(T, Exp1Par)

           Exp1_TMP = exp(Exp1Par(1) * T^2 + Exp1Par(2) * T + Exp1Par(3));

        end

        %% Evaluate the Kim and Lee for plotting purposes

        function KL_TMP = KL_Plot(T, KLPar)

            KL_TMP = 1 - (KLPar(1) * exp(1 + ((KLPar(2) - T) / KLPar(3)) - ...
                          exp((KLPar(2) - T) / KLPar(3))));

        if KL_TMP < 0

            KL_TMP = 0;
    
        else

            KL_TMP = abs(1 - (KLPar(1) * exp(1 + ((KLPar(2) - T) / KLPar(3)) - ...
                          exp((KLPar(2) - T) / KLPar(3)))));
        end

        end


        %% Evaluate the Gaussian-like for plotting purposes

        function Fert_TMP = Fert_Plot(T, FertPar)

            Fert_TMP = FertPar(1) * ( ((FertPar(2) + 1) / (pi * ...
                       (FertPar(3) ^ (2 * FertPar(2) + 2)))) *  ...
                       ((FertPar(3) ^ 2) - ( ((T - FertPar(5)) ^ 2) ...
                       + (FertPar(4) ^ 2)) ) ^ FertPar(2) );

        end

    end

end





