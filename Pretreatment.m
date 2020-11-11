%% MATLAB Thermodynamics of Pretreatment Step of Biodiesel Production [REA]
%
% Members:
%
% Owen Wilborn
% Rebecca Soares
%
% Project start date 10/22/2020
%
% This code simulations the conversion of different free fatty acids [FFA] into
% free fatty mythel esters [FFME] by analyzing the thermodynmaics of the
% the rxn [FFA + Methonal -> FFME + H2O]

%%  Clear Cache

clear
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Raw Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input all data in standard SI units of M,kg,s,N,J,K,Mol

%% Standard Variables

R = 8.3144626; %Ideal gas constant (J/(K*Mol))
T0 = 298.15; %K, standard temperature

%% Class Molecule

%% Methanol (CH3OH)

%Create molecule
methonal = Molecule;
methonal.name = "methonal";

%Standard Formation Energies
methonal.Enthalpy_Formation = -239100; %J/mol, liquid
methonal.Gibbs_Formation = -166900; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 175.47K to 503.15K
methonal.A_Heat_Capacity = 256040;
methonal.B_Heat_Capacity = -2741.4;
methonal.C_Heat_Capacity = 14.777;
methonal.D_Heat_Capacity = -0.035078;
methonal.E_Heat_Capacity = 3.2719*10^-5;

%% Water (H2O)

%Create molecule
water = Molecule;
water.name = "water";

%Standard Formation Energies
water.Enthalpy_Formation = -285830; %J/mol, liquid
water.Gibbs_Formation = -237129; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 273.16K to 533.15K
water.A_Heat_Capacity = 276370;
water.B_Heat_Capacity = -2090.1;
water.C_Heat_Capacity = 8.125;
water.D_Heat_Capacity  = -0.014116;
water.E_Heat_Capacity = 9.3701*10^-6;

%% Oleic Acid (cis-9-octadecenoic acid)

%Create molecule
oleic_acid = Molecule;
oleic_acid.name = "oleic acid";

%Standard Formation Energies
oleic_acid.Enthalpy_Formation = -802491; %J/mol, liquid
oleic_acid.Gibbs_Formation = -261100; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 286.53K to 550K
oleic_acid.A_Heat_Capacity = 459000;
oleic_acid.B_Heat_Capacity = -866;
oleic_acid.C_Heat_Capacity = 3.74;
oleic_acid.D_Heat_Capacity = 0;
oleic_acid.E_Heat_Capacity = 0;


%% Methyl Oleate (methyl (Z)-octadec-9-enoate)

%Create molecule
methyl_oleate = Molecule;
methyl_oleate.name = "methyl oleate";

%Standard Formation Energies
methyl_oleate.Enthalpy_Formation = -734500; %J/mol, liquid
methyl_oleate.Gibbs_Formation = -170900; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 293.53K to 617K
methyl_oleate.A_Heat_Capacity = 324000;
methyl_oleate.B_Heat_Capacity = 928;
methyl_oleate.C_Heat_Capacity = 0;
methyl_oleate.D_Heat_Capacity = 0;
methyl_oleate.E_Heat_Capacity = 0;

%% Linoleic Acid ((9Z,12Z)-octadeca-9,12-dienoic acid)

%Create molecule
linoleic_acid = Molecule;
linoleic_acid.name = "linoleic acid";

%Standard Formation Energies
linoleic_acid.Enthalpy_Formation = -674042; %J/mol, liquid
linoleic_acid.Gibbs_Formation = -168200; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 387.5K to 581K
linoleic_acid.A_Heat_Capacity = 122700;
linoleic_acid.B_Heat_Capacity = 1214.5;
linoleic_acid.C_Heat_Capacity = 0;
linoleic_acid.D_Heat_Capacity = 0;
linoleic_acid.E_Heat_Capacity = 0;

%% Methyl Linoleate (methyl (9Z,12Z)-octadeca-9,12-dienoate)

%Create molecule
methyl_linoleate = Molecule;
methyl_linoleate.name = "methyl linoleate";

%Standard Formation Energies
methyl_linoleate.Enthalpy_Formation = -661700; %J/mol, liquid
methyl_linoleate.Gibbs_Formation = -101800; %J/mol, liquid

%Heat capacity of liquid
% Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4 for T 233.9K to 690K
methyl_linoleate.A_Heat_Capacity = 552150;
methyl_linoleate.B_Heat_Capacity = -707.45;
methyl_linoleate.C_Heat_Capacity = 2.9016;
methyl_linoleate.D_Heat_Capacity = 0;
methyl_linoleate.E_Heat_Capacity = 0;

%% %%%%%%%%%%%%%%% 1. Pretreatment Equilibruim Conversion Modeling %%%%%%%%%%%%%%%%%%%%%%%%

%% Define each reaction

% Oleic acid -> Mythel Oleate
rxn_oleic = Reaction;

rxn_oleic.Reactant(1) = methonal;
rxn_oleic.Reactant(2) = oleic_acid;
rxn_oleic.Product(1)= water;
rxn_oleic.Product(2) = methyl_oleate;

setRxn(rxn_oleic);

% Linoleic acid -> Mythel linolate
rxn_linoleic = Reaction;

rxn_linoleic.Reactant(1) = methonal;
rxn_linoleic.Reactant(2) = linoleic_acid;
rxn_linoleic.Product(1)= water;
rxn_linoleic.Product(2) = methyl_linoleate;

setRxn(rxn_linoleic);

%% Equilbrium conversion at 298.15K
K0_Oleic = exp(-rxn_oleic.Delta_Gibbs/(R*T0));
K0_Linoleic = exp(-rxn_linoleic.Delta_Gibbs/(R*T0));

%% Equlibrum conversion dependence on temperature

%Set counters
Told = 298.15;
CP_Heat_Affect_Oleic = 0;
CP_Heat_Affect_Linoleic = 0;
i=1;
syms temp

%For loop to iterate through multiple temperatures
for T=298.15:1:398.15
    
    %K1
    K1_Oleic = exp(rxn_oleic.Delta_Enthalpy/(R*T0)*(1-T0/T));
    K1_Linoleic = exp(rxn_linoleic.Delta_Enthalpy/(R*T0)*(1-T0/T));  

    %K2
    CP_Oleic = get_CP_Function(rxn_oleic);
    CP_Heat_Addition_Oleic = (-1/(T*R)*(int((CP_Oleic), temp, Told, T))+ 1/R*(int((CP_Oleic/temp), temp, Told, T)));
    CP_Heat_Affect_Oleic = CP_Heat_Affect_Oleic + CP_Heat_Addition_Oleic;
    K2_Oleic = exp(CP_Heat_Affect_Oleic);
    
    CP_Linoleic = get_CP_Function(rxn_linoleic);
    CP_Heat_Addition_Linoleic = (-1/(T*R)*(int((CP_Linoleic), temp, Told, T))+ 1/R*(int((CP_Linoleic/temp), temp, Told, T)));
    CP_Heat_Affect_Linoleic = CP_Heat_Affect_Linoleic + CP_Heat_Addition_Linoleic;
    K2_Linoleic = exp(CP_Heat_Affect_Linoleic);
    
    %Solve for equilibrium conversion K and save data
    K_oleic(i) = K0_Oleic*K1_Oleic*K2_Oleic;
    K_linoleic(i) = K0_Linoleic* K1_Linoleic * K2_Linoleic;
    Temp(i) = T;
    
    %Increase Step size and counters
    i=i+1;
    Told = T;
end

%Change data to decimal form

K_oleic = double(K_oleic);
K_linoleic = double(K_linoleic);

%Plot conversion K vs. Temp.

figure(1)
plot(Temp, (K_oleic))
title('Effect of Temperature on Conversion of Oleic Acid to Methyl Oleate')
ylabel('Equilibrium Constant K') 
ax = gca;
ax.YAxis.Exponent = 0;
ytickformat('%.1e')
xlabel('Temperature (K)')

figure(2)
plot(Temp, (K_linoleic))
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
        Error_Oleic = K_oleic(i) - ((N_int_Water+v_Water*extent_rxn)*(N_int_Methyl_Oleate+v_Methyl_Oleate*extent_rxn)/((N_int_Oleic_Acid+v_Oleic_Acid*extent_rxn)*(N_int_Methonal+v_Methonal*extent_rxn)));
        if abs(Error_Oleic) < abs(errorMax_Oleic)
            extent_rxn_Oleic(i) = extent_rxn;
            errorMax_Oleic = Error_Oleic;
        end
        
        Error_Linoleic = K_linoleic(i) - ((N_int_Water+v_Water*extent_rxn)*(N_int_Methyl_Linoleate+v_Methyl_Linoleate*extent_rxn)/((N_int_Linoleic_Acid+v_Linoleic_Acid*extent_rxn)*(N_int_Methonal+v_Methonal*extent_rxn)));
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




