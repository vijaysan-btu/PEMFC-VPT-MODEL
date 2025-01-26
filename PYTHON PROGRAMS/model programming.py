#Library import
import math
import numpy

#Variables
T=float(input("Enter Operating temperature : "))
Ph2=float(input("Enter Operating pressure of Hydrogen : "))
PO2=float(input("Enter Operating pressure of Oxygen : "))
A=float(input("Enter the Active area : "))
i = float(input("Enter the starting current : "))
lamb = int(inpit("Enter the value of Lambda(14-23) : ")) #Mostly consider 23
I=float(input("Enter the Membrane thickness : "))
N = int(input("Enter the number of cells : "))
float Cell_Voltage,Reference_Volatge,Nernst_Voltage,Stack_Voltage

#Nernst Volatge Calculation
Reference_Volatage=1.229
Nernst_Volatge=Reference_Volatge - ((8.5*e-4)*(T-298.15)) + ((4.308*e-5)*T*((math.log(Ph2,10))+(math.log(PO2)))

#Losses Calculation

#Activation losses Calculation
float yeta_1=-0.948,yeta_2,yeta_3=7.6*e-5,yeta_4=-1.93*e-4
float Ch2 = Ph2/(1.09*e6*(math.exp(77/T)))
float Co2 = PO2/(5.08*e6*(math.exp(-498/T)))
yeta_2=0.00286 + 0.002*(math.log(A,10)) + (4.3*e-5)*(math.log(Ch2,10))
float Activation-loss = yeta_1 + (yeta_2 * T) + (yeta_3 * T * (math.log(Co2))) + (yeta_4 * T * (math.log(i)))

#Ohmic losses Calculation
int R_electronic = 0 #Assuming losses due to Electronic components to be negligible
float c-a = i/A
float ep = (4.18*((T-303)/T))
float Rho = (181.6(1 + (0.03*c-a) + (0.062*math.pow(T/303,2)*math.pow(c-a,2.5))))/(lamb - 0.634 - (3*c-a*math.exp(ep)))
float R-Proton = (Rho * I)/A 
float Ohmic-Loss = i*(R_Electronic + R-proton)


#Concentration loss Calculation
float B = 0.015
float J-max = 1.5
float J = c-a
float Concentration-loss = -B * math.log((1-(J/J-max)),10)

#Total loss calculation
float Loss = Activation-loss + Ohmic-loss + Concentration-loss

#Voltage calculation
Cell_Voltage = Nernst_Volatge - Loss
Stack_Voltage = N * Cell_Voltage

float P-cell = Cell_Volatge * i #Powe Calculation
float P-Stack = N * P-cell