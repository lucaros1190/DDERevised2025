# Temperature-dependent delays - To feed the DDE

# Import list

from Parameters import *


# Temperature reader

def TempFunction(time, TempArray):

    if time >= len(TempArray):
        Temp = TempArray[0]
    else:
        Temp = TempArray[time]

    return Temp


# Delay function 1 - 2nd order rational polynomial

def Delayfun_RatPoly(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]
    d = Parameters[3]

    Delay = (a * T**2 + b * T + c) / (T + d)

    return Delay


# Delay function 2 - 2nd order polynomial (parabola)

def Delayfun_Parab(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = a * T**2 + b * T + c

    return Delay


# Delay function 3 - Exponential

def Delayfun_Expo(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = a + b * np.exp(c * T)

    return Delay


# Temperature-dependent delays - To feed the DDE

def Delays(t, TempArray, LagPar_Egg, LagPar_L1, LagPar_L2, LagPar_L3, LagPar_P, LagPar_AF2):
    
    # Manage time and temperatures

    time = int(t)

    DayTemp = TempFunction(time, TempArray)

        # Egg lag
    tau_E = Delayfun_RatPoly(DayTemp, LagPar_Egg)

        # L1 lag
    tau_L1 = Delayfun_Parab(DayTemp, LagPar_L1)

        # L2 lag
    tau_L2 = Delayfun_Parab(DayTemp, LagPar_L2)

        # L3 lag
    tau_L3 = Delayfun_Parab(DayTemp, LagPar_L3)

        # P_lag
    tau_P = Delayfun_Expo(DayTemp, LagPar_P)

        # Male
    tau_AM = 0
    
        # Female 1
    tau_AF1 = 0

        # Female 2 lag (PreOvi)
    tau_AF2 = Delayfun_Parab(DayTemp, LagPar_AF2)

    Del = [tau_E, tau_L1, tau_L2, tau_L3, tau_P, tau_AM, tau_AF1, tau_AF2]

    return np.array(Del) / 1.8
    

# Definition of Logan rate function
    
def Log(Par, T):

    G = Par[0] * (np.exp(Par[1] * T) - np.exp(Par[1] * Par[2] - ((Par[2] - T) / Par[3])))
    
    if G > 0:
        G = Par[0] * (np.exp(Par[1] * T) - np.exp(Par[1] * Par[2] - ((Par[2] - T) / Par[3])))
    else:
        G = 0
        
    return G


# Definition of Sharpe and De Michele rate function

def SDM(Par, Temp):

    T = Temp + 273.15 # Temperature in Kelvin for this function!
    
    if T <= 0:
        T = 1
    else:
        T = Temp

    G = (T * (np.exp(Par[0] - (Par[1] / T)))) / ((1 + np.exp(Par[2] - (Par[3] / T)) + np.exp(Par[4] - (Par[5] / T))))

    if G > 0:
        G = (T * (np.exp(Par[0] - (Par[1] / T)))) / ((1 + np.exp(Par[2] - (Par[3] / T)) + np.exp(Par[4] - (Par[5] / T))))
    else:
        G = 0
        
    return G


# Definition of Lactin-1 rate function

def LacOne(Par, T):

    G = np.exp(Par[0] * T) - np.exp(Par[0] * Par[1] - ( (Par[1] - T) / Par[2] ) )

    if G > 0:
        G = np.exp(Par[0] * T) - np.exp(Par[0] * Par[1] - ( (Par[1] - T) / Par[2] ) )
    else:
        G = 0
        
    return G
    
    
# Definition of Briere rate function

def Bri(Par, T):

    G = Par[0] * T * (T - Par[1]) * np.pow(Par[2] - T, (1/Par[3]))

    if G > 0:
        G = Par[0] * T * (T - Par[1]) * np.pow(Par[2] - T, (1/Par[3]))
    else:
        G = 0
        
    return G
    

# Definition of the fertility rate function

def FertFunc(Par, T):

    FP = Par[0] * ( ((Par[1] + 1) / (np.pi * (Par[2] ** (2 * Par[1] + 2)))) * ((Par[2] ** 2) - ( ((T - Par[4]) ** 2) + (Par[3] ** 2)) ) ** Par[1] )
    
    return 0# FP


# Definition of the Mortality rate function - Exponential with parabola

def MortFunc(Par, T):

    MP = np.exp(Par[0] * T**2 + Par[1] * T + Par[2])

    if MP >= 0:
        MP = np.exp(Par[0] * T ** 2 + Par[1] * T + Par[2])
    else:
        MP = 0

    return MP


# Compute the normalisation coefficients for the development rates – Including adult survival

def NormCoefficients(DevPar_Egg, DevPar_L1, DevPar_L2, DevPar_L3,
                     DevPar_P, DevPar_AF, DevPar_AM):
                    
    # Array to store solutions
    
    DevRate_Egg = []
    DevRate_L1 = []
    DevRate_L2 = []
    DevRate_L3 = []
    DevRate_P = []
    DevRate_AF = []
    DevRate_AM = []
    
    # Temperature series to evaluate the functions in their domain ox existance
    
    Temp = np.linspace(1, 50, 1000)
    
    # Compute development rates in their domain of existance
    
    DevRate_Egg = np.array([SDM(DevPar_Egg, T) for T in Temp])
    DevRate_L1 = np.array([Bri(DevPar_L1, T) for T in Temp])
    DevRate_L2 = np.array([Log(DevPar_L2, T) for T in Temp])
    DevRate_L3 = np.array([LacOne(DevPar_L3, T) for T in Temp])
    DevRate_P = np.array([Bri(DevPar_P, T) for T in Temp])
    DevRate_AF = np.array([LacOne(DevPar_AF, T) for T in Temp])
    DevRate_AM = np.array([LacOne(DevPar_AM, T) for T in Temp])
    
    # Compute maximum values
    
    Egg_Norm = np.max(DevRate_Egg)
    L1_Norm = np.max(DevRate_L1)
    L2_Norm = np.max(DevRate_L2)
    L3_Norm = np.max(DevRate_L3)
    P_Norm = np.max(DevRate_P)
    AF_Norm = np.max(DevRate_AF)
    AM_Norm = np.max(DevRate_AM)
    
    return [Egg_Norm, L1_Norm, L2_Norm, L3_Norm, P_Norm, AF_Norm, AM_Norm]
    

# Compute Egg daily development rate normalised

def EggNormDevRate(Par, T, NormCoeff):

    EggNormDevRate_Value = SDM(Par, T) #/ NormCoeff
    
    return EggNormDevRate_Value * 2


# Compute L1 daily development rate normalised

def L1NormDevRate(Par, T, NormCoeff):

    L1NormDevRate_Value = Bri(Par, T) #/ NormCoeff
    
    return L1NormDevRate_Value  * 2


# Compute L2 daily development rate normalised

def L2NormDevRate(Par, T, NormCoeff):

    L2NormDevRate_Value = Log(Par, T) #/ NormCoeff

    return L2NormDevRate_Value * 2


# Compute L3 daily development rate normalised

def L3NormDevRate(Par, T, NormCoeff):

    L3NormDevRate_Value = LacOne(Par, T) #/ NormCoeff

    return L3NormDevRate_Value * 2
    
    
# Compute Pupa daily development rate normalised

def PNormDevRate(Par, T, NormCoeff):

    PNormDevRate_Value = Bri(Par, T) #/ NormCoeff
    
    return PNormDevRate_Value * 2


# Compute Adult Female daily development rate normalised

def AFNormDevRate(Par, T, NormCoeff):

    AFNormDevRate_Value = LacOne(Par, T) #/ NormCoeff
    
    return AFNormDevRate_Value
    

# Compute Adult Male daily development rate normalised

def AMNormDevRate(Par, T, NormCoeff):

    AMNormDevRate_Value = LacOne(Par, T) #/ NormCoeff
    
    return AMNormDevRate_Value
    

# DDE system REVISED - To feed the DDE solver!

def DDE_FuncNormRevised(Y, t, FertPar, MortPar_Egg, MortPar_L1, MortPar_L2,
                        MortPar_L3, MortPar_P, DevPar_Egg, DevPar_L1, DevPar_L2,
                        DevPar_L3, DevPar_P, DevPar_AF, DevPar_AM,
                        LagPar_Egg, LagPar_L1, LagPar_L2, LagPar_L3,
                        LagPar_P, LagPar_AF2, SR, TempArray, NormCoeff_Egg,
                        NormCoeff_L1, NormCoeff_L2, NormCoeff_L3, NormCoeff_P, 
                        NormCoeff_AF, NormCoeff_AM):
                   
    # Calculate the daily temperature from the 'TemperatureInput.xlsx' file in 'DS_Parameters.py'

    temp = TempFunction(int(t), TempArray)

    # Approximation of the lags - Needed by ddeint in Y(t - Z)

    # Z[0] = Egg lag
    # Z[1] = L1 lag
    # Z[2] = L2 lag
    # Z[3] = L3 lag
    # Z[4] = P lag
    # Z[5] = AM
    # Z[6] = AF1
    # Z[7] = AF2 lag (Preoviposition)

        # Assigning to variables for not becoming crazy on the following part of the code!

    Z = Delays(t, TempArray, LagPar_Egg, LagPar_L1, LagPar_L2, LagPar_L3, LagPar_P, LagPar_AF2)

    tau_E = Z[0]
    tau_L1 = Z[1]
    tau_L2 = Z[2]
    tau_L3 = Z[3]
    tau_P = Z[4]
    tau_AM = Z[5]
    tau_AF1 = Z[6]
    tau_AF2 = Z[7]

        # Checking the values for the initial history

    Egg_tau = Y(t - tau_E)[0]
    L1_tau = Y(t - tau_L1)[2]
    L2_tau = Y(t - tau_L2)[4]
    L3_tau = Y(t - tau_L3)[6]
    P_tau = Y(t - tau_P)[8]
    AM_tau = Y(t - tau_AM)[10]
    AF1_tau = Y(t - tau_AF1)[11]
    AF2_tau = Y(t - tau_AF2)[12]

    # DDE system

    # Y(t)[0] = Egg without lag
    # Y(t)[1] = Transient eggs
    # Y(t)[2] = L1 without lag
    # Y(t)[3] = Transient L1
    # Y(t)[4] = L2 without lag
    # Y(t)[5] = Transient L2
    # Y(t)[6] = L3 without lag
    # Y(t)[7] = Transient L3
    # Y(t)[8] = P without lag
    # Y(t)[9] = Transient P
    # Y(t)[10] = Male without lag
    # Y(t)[11] = AF1 without lag
    # Y(t)[12] = AF2 without lag

    Egg = Y(t)[0]
    Egg_Tr = Y(t)[1]
    L1 = Y(t)[2]
    L1_Tr = Y(t)[3]
    L2 = Y(t)[4]
    L2_Tr = Y(t)[5]
    L3 = Y(t)[6]
    L3_Tr = Y(t)[7]
    P = Y(t)[8]
    P_Tr = Y(t)[9]
    AM = Y(t)[10]
    AF1 = Y(t)[11]
    AF2 = Y(t)[12]
    
    # Pre-computing mortalities
    
    M_E = MortFunc(MortPar_Egg, temp) * SDM(DevPar_Egg, temp)
    M_Etr = MortFunc(MortPar_Egg, temp) * SDM(DevPar_Egg, temp)
    
    M_L1 = MortFunc(MortPar_L1, temp) * Bri(DevPar_L1, temp)
    M_L1tr = MortFunc(MortPar_L1, temp) * Bri(DevPar_L1, temp)
    
    M_L2 = MortFunc(MortPar_L2, temp) * Log(DevPar_L2, temp)
    M_L2tr = MortFunc(MortPar_L2, temp) * Log(DevPar_L2, temp)
    
    M_L3 = MortFunc(MortPar_L3, temp) * LacOne(DevPar_P, temp)
    M_L3tr = MortFunc(MortPar_L3, temp) * LacOne(DevPar_P, temp)
    
    M_P = MortPar_P
    M_Ptr = MortPar_P
    
    
        # Egg

    dEgg = LacOne(DevPar_AF, temp) * FertFunc(FertPar, temp) * AF2_tau - EggNormDevRate(DevPar_Egg, temp, NormCoeff_Egg) * Egg - M_E * Egg

        # Transient Egg

    dEgg_Tr = EggNormDevRate(DevPar_Egg, temp, NormCoeff_Egg) * Egg - np.exp(- M_Etr * tau_E) * EggNormDevRate(DevPar_Egg, temp, NormCoeff_Egg) * Egg_tau - M_Etr * Egg_Tr

        # L1
    
    dL1 = np.exp(- M_Etr * tau_E) * EggNormDevRate(DevPar_Egg, temp, NormCoeff_Egg) * Egg_tau - L1NormDevRate(DevPar_L1, temp, NormCoeff_L1) * L1 - M_L1 * L1

        # Transient L1
    
    dL1_Tr = L1NormDevRate(DevPar_L1, temp, NormCoeff_L1) * L1 - np.exp(- M_L1tr * tau_L1) * L1NormDevRate(DevPar_L1, temp, NormCoeff_L1) * L1_tau - M_L1tr * L1_Tr

        # L2
    
    dL2 = np.exp(- M_L1tr * tau_L1) * L1NormDevRate(DevPar_L1, temp, NormCoeff_L1) * L1_tau - L2NormDevRate(DevPar_L2, temp, NormCoeff_L2) * L2 - M_L2 * L2

        # Transient L2
    
    dL2_Tr = L2NormDevRate(DevPar_L2, temp, NormCoeff_L2) * L2 - np.exp(- M_L2tr * tau_L2) * L2NormDevRate(DevPar_L2, temp, NormCoeff_L2) * L2_tau - M_L2tr * L2_Tr

        # L3
    
    dL3 = np.exp(- M_L2tr * tau_L2) * L2NormDevRate(DevPar_L2, temp, NormCoeff_L2) * L2_tau - L3NormDevRate(DevPar_L3, temp, NormCoeff_L3) * L3 - M_L3 * L3

        # Transient L3
    
    dL3_Tr = L3NormDevRate(DevPar_L3, temp, NormCoeff_L3) * L3 - np.exp(- M_L3tr * tau_L3) * L3NormDevRate(DevPar_L3, temp, NormCoeff_L3) * L3_tau - M_L3tr * L3_Tr

        # P
    
    dP = np.exp(- M_L3tr * tau_L3) * L3NormDevRate(DevPar_L3, temp, NormCoeff_L3) * L3_tau - PNormDevRate(DevPar_P, temp, NormCoeff_P) * P - M_P * P

        # Transient P
    
    dP_Tr =  PNormDevRate(DevPar_P, temp, NormCoeff_P) * P - np.exp(- M_Ptr * tau_P) * PNormDevRate(DevPar_P, temp, NormCoeff_P) * P_tau - M_Ptr * P_Tr

        # AM
    
    dAM = np.exp(- M_Ptr * tau_P) * (1 - SR) * PNormDevRate(DevPar_P, temp, NormCoeff_P) * P_tau - (AMNormDevRate(DevPar_AM, temp, NormCoeff_AM)) * AM

        # Female-1 (Non mated)
    
    dAF1 = np.exp(- M_Ptr * tau_P) * SR * PNormDevRate(DevPar_P, temp, NormCoeff_P) * P_tau - AF1

        # Female-2 (Mated)
    
    dAF2 = AF1 - (AFNormDevRate(DevPar_AF, temp, NormCoeff_AF)) * AF1 - (AFNormDevRate(DevPar_AF, temp, NormCoeff_AF)) * AF2
    
    dydt_partial = [dEgg, dEgg_Tr, dL1, dL1_Tr, dL2, dL2_Tr, dL3, dL3_Tr, dP, dP_Tr, dAM, dAF1, dAF2]

    return np.array(dydt_partial)
