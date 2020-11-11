%% Class Definition for Molecule
%
% Molecule class for each specific molecule used. Each molecule has an
% enthalpy and gibbs energy of formations and A,B,C,D,E heat capacity
% coefficients.
%
classdef Molecule 
    properties
        % Name
        name = string();%Stantard name of the molecule
        
        % Thermodynamic Properties
        
        Enthalpy_Formation {mustBeNumeric} %Standard Enthalpy of Formation @ 298 K
        Gibbs_Formation {mustBeNumeric} %Standard Gibbs Energy of Formation @ 298 K
        
        % Heat Capacity Liquid
        % Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4
        A_Heat_Capacity {mustBeNumeric}
        B_Heat_Capacity {mustBeNumeric}
        C_Heat_Capacity {mustBeNumeric}
        D_Heat_Capacity {mustBeNumeric}
        E_Heat_Capacity {mustBeNumeric}
        
    end
end