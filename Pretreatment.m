%% MATLAB Thermodynamics of Pretreatment Step of Biodiesel Production [REA]

%%  Clear Cache
clear all
close all
clc

%% Raw Data
%% Standard Variables
R = 8.3144626; %Ideal gas constant (J/(K*Mol))
T0 = 298.15; %K, standard temperature

%% Methanol

%Standard Formation Energies
Hf_Methonal = -239100; %J/mol, liquid
Gf_Methonal = -166900; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 175.47K to 503.15K
A_Methonal = 256040;
B_Methonal = -2741.4;
C_Methonal = 14.777;
D_Methonal = -0.035078;
E_Methonal = 3.2719*10^-5;

%% Water

%Standard Formation Energies
Hf_Water = -285830; %J/mol, liquid
Gf_Water = -237129; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 273.16K to 533.15K
A_Water = 276370;
B_Water = -2090.1;
C_Water = 8.125;
D_Water  = -0.014116;
E_Water = 9.3701*10^-6;

%% Oleic Acid (cis-9-octadecenoic acid)

%Standard Formation Energies
Hf_Oleic_Acid = -802491; %J/mol, liquid
Gf_Oleic_Acid = -261100; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 286.53K to 550K
A_Oleic_Acid = 459000;
B_Oleic_Acid = -866;
C_Oleic_Acid = 3.74;
D_Oleic_Acid = 0;
E_Oleic_Acid = 0;


%% Methyl Oleate (methyl (Z)-octadec-9-enoate)

%Standard Formation Energies
Hf_Methyl_Oleate = -734500; %J/mol, liquid
Gf_Methyl_Oleate = -170900; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 293.53K to 617K
A_Methyl_Oleate = 324000;
B_Methyl_Oleate = 928;
C_Methyl_Oleate = 0;
D_Methyl_Oleate = 0;
E_Methyl_Oleate = 0;

%% Linoleic Acid ((9Z,12Z)-octadeca-9,12-dienoic acid)

%Standard Formation Energies
Hf_Linoleic_Acid = -674042; %J/mol, liquid
Gf_Linoleic_Acid = -168200; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 387.5K to 581K
A_Linoleic_Acid = 122700;
B_Linoleic_Acid = 1214.5;
C_Linoleic_Acid = 0;
D_Linoleic_Acid = 0;
E_Linoleic_Acid = 0;

%% Methyl Linoleate (methyl (9Z,12Z)-octadeca-9,12-dienoate)

%Standard Formation Energies
Hf_Methyl_Linoleate = -661700; %J/mol, liquid
Gf_Methyl_Linoleate = -101800; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 233.9K to 690K
A_Methyl_Linoleate = 552150;
B_Methyl_Linoleate = -707.45;
C_Methyl_Linoleate = 2.9016;
D_Methyl_Linoleate = 0;
E_Methyl_Linoleate = 0;

%% Reaction Modeling

%% Rxn 1 Pretreatment Equilibruim Conversion

% Formation energies

% [Methanol + Oleic Acid -> Methyl Olate + Water]
delta_G_Oleic = (Gf_Methyl_Oleate + Gf_Water) - (Gf_Oleic_Acid + Gf_Methonal); %(J/mol)
delta_H_Oleic = (Hf_Methyl_Oleate + Hf_Water) - (Hf_Oleic_Acid + Hf_Methonal); %(J/mol)

% [Methanol + Linoleic Acid -> Linoleate Olate + Water]
delta_G_Linoleic = (Gf_Methyl_Linoleate + Gf_Water) - (Gf_Linoleic_Acid + Gf_Methonal); %(J/mol)
delta_H_Linoleic = (Hf_Methyl_Linoleate + Hf_Water) - (Hf_Linoleic_Acid + Hf_Methonal); %(J/mol)

%Equilbrium conversion at 298.15K

%K0

K0_Oleic = exp(-delta_G_Oleic/(R*T0));
K0_Linoleic = exp(-delta_G_Linoleic/(R*T0));

%Conversion at varying temperature

A_Oleic = (A_Methyl_Oleate + A_Water) - (A_Oleic_Acid + A_Methonal);
B_Oleic = (B_Methyl_Oleate + B_Water) - (B_Oleic_Acid + B_Methonal);
C_Oleic = (C_Methyl_Oleate + C_Water) - (C_Oleic_Acid + C_Methonal);
D_Oleic = (D_Methyl_Oleate + D_Water) - (D_Oleic_Acid + D_Methonal);
E_Oleic = (E_Methyl_Oleate + E_Water) - (E_Oleic_Acid + E_Methonal);
syms temp
CP_Oleic = @(temp) ((A_Oleic + B_Oleic*temp + C_Oleic*temp^2 + D_Oleic*temp^3 + E_Oleic*temp^4)/1000); %(J/(mol*K)

A_Linoleic = (A_Methyl_Linoleate + A_Water) - (A_Linoleic_Acid + A_Methonal);
B_Linoleic = (B_Methyl_Linoleate + B_Water) - (B_Linoleic_Acid + B_Methonal);
C_Linoleic = (C_Methyl_Linoleate + C_Water) - (C_Linoleic_Acid + C_Methonal);
D_Linoleic = (D_Methyl_Linoleate + D_Water) - (D_Linoleic_Acid + D_Methonal);
E_Linoleic = (E_Methyl_Linoleate + E_Water) - (E_Linoleic_Acid + E_Methonal);
syms temp
CP_Linoleic = @(temp) ((A_Linoleic + B_Linoleic*temp + C_Linoleic*temp^2 + D_Linoleic*temp^3 + E_Linoleic*temp^4)/1000); %(J/(mol*K)

%Set counters
Told = 298.15;
CP_Heat_Affect_Oleic = 0;
CP_Heat_Affect_Linoleic = 0;
i=1;

%For loop to iterate through multiple temperatures
for T=298.15:1:398.15
    
    %K1
    K1_Oleic = exp(delta_H_Oleic/(R*T0)*(1-T0/T));
    K1_Linoleic = exp(delta_H_Linoleic/(R*T0)*(1-T0/T));  

    %K2
    CP_Heat_Addition_Oleic = (-1/(T*R)*(int((CP_Oleic), temp, Told, T))+ 1/R*(int((CP_Oleic/temp), temp, Told, T)));
    CP_Heat_Affect_Oleic = CP_Heat_Affect_Oleic + CP_Heat_Addition_Oleic;
    K2_Oleic = exp(CP_Heat_Affect_Oleic);
    
    CP_Heat_Addition_Linoleic = (-1/(T*R)*(int((CP_Linoleic), temp, Told, T))+ 1/R*(int((CP_Linoleic/temp), temp, Told, T)));
    CP_Heat_Affect_Linoleic = CP_Heat_Affect_Linoleic + CP_Heat_Addition_Linoleic;
    K2_Linoleic = exp(CP_Heat_Affect_Linoleic);
    
    %Solve for equilibrium conversion K and save data
    K_Oleic(i) = K0_Oleic*K1_Oleic*K2_Oleic;
    K_Linoleic(i) = K0_Linoleic* K1_Linoleic * K2_Linoleic;
    Temp(i) = T;
    
    %Increase Step size and counters
    i=i+1;
    Told = T;
end

K_Oleic = double(K_Oleic);
K_Linoleic = double(K_Linoleic);

%Plot conversion K vs. Temp.

figure(1)
plot(Temp, (K_Oleic))
title('Effect of Temperature on Conversion of Oleic Acid to Methyl Oleate')
ylabel('Equilibrium Constant K') 
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.1e')
xlabel('Temperature (K)')

figure(2)
plot(Temp, (K_Linoleic))
title('Effect of Temperature on Conversion of Linoleic Acid to Methyl Linoleate')
ylabel('Equilibrium Constant K') 
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.1e')
xlabel('Temperature (K)')

%Calculate Conversion at each temperature

N_int_Methonal = 20;         v_Methonal = -1;
N_int_Oleic_Acid = 1;         v_Oleic_Acid = -1;
N_int_Methyl_Oleate = 0;     v_Methyl_Oleate = 1;
N_int_Linoleic_Acid = 1;              v_Linoleic_Acid = -1;
N_int_Methyl_Linoleate = 0;     v_Methyl_Linoleate = 1;
N_int_Water = 0;                v_Water = 1;


%Find extent of rxn
for i=1:1:101
   errorMax_Oleic = 100;
   errorMax_Linoleic = 100;
   
   for extent_rxn = 0:0.0001:1
        Error_Oleic = K_Oleic(i) - ((N_int_Water+v_Water*extent_rxn)*(N_int_Methyl_Oleate+v_Methyl_Oleate*extent_rxn)/((N_int_Oleic_Acid+v_Oleic_Acid*extent_rxn)*(N_int_Methonal+v_Methonal*extent_rxn)));
        if abs(Error_Oleic) < abs(errorMax_Oleic)
            extent_rxn_Oleic(i) = extent_rxn;
            errorMax_Oleic = Error_Oleic;
        end
        
        Error_Linoleic = K_Linoleic(i) - ((N_int_Water+v_Water*extent_rxn)*(N_int_Methyl_Linoleate+v_Methyl_Linoleate*extent_rxn)/((N_int_Linoleic_Acid+v_Linoleic_Acid*extent_rxn)*(N_int_Methonal+v_Methonal*extent_rxn)));
        if abs(Error_Linoleic) < abs(errorMax_Linoleic)
            extent_rxn_Linoleic(i) = extent_rxn;
            errorMax_Linoleic = Error_Linoleic;
        end
    end
end

%Plot extent of reaction for different temperatures
figure(3)
hold on
plot(Temp, extent_rxn_Oleic)
plot(Temp, extent_rxn_Linoleic)
title('Effect of Temperature on Extent of Reaction')
ylabel('Extent of Reaction') 
xlabel('Temperature (K)')
legend('Oleic Acid', 'Linoleic Acid')
hold off




