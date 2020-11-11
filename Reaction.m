%% Class Definition for Reaction
%
%Class for a reaction when create must add the reactants and products. Once
%calls the function set reaction the class variables for deltaG, deltaH and
%A,B,C,D,E for heat cap. are calculated and set
%
classdef Reaction < handle
    %% Properties
    properties
        % Vector of molecules on the products side of the equations
        Product = Molecule();
        
        % Vector of molecules on the reactant side of the equations
        Reactant  = Molecule();
        
        % Thermodynamic Properties of the reaction
        
        Delta_Enthalpy {mustBeNumeric} %Standard Enthalpy of Formation @ 298 K
        Delta_Gibbs {mustBeNumeric} %Standard Gibbs Energy of Formation @ 298 K
        
        % Heat Capacity Liquid
        % Cp (J/(kmol*K)) = A + BT + CT^2 + DT^3 + ET^4
        A_Heat_Capacity {mustBeNumeric}
        B_Heat_Capacity {mustBeNumeric}
        C_Heat_Capacity {mustBeNumeric}
        D_Heat_Capacity {mustBeNumeric}
        E_Heat_Capacity {mustBeNumeric}
    end
    
    %% Public Methods
    methods (Access=public)
       function deltaG = get_deltaG(obj)
           deltaG = obj.Delta_Gibbs;
       end
       function deltaH = get_deltaH(obj)
           deltaH = obj.Delta_Enthalpy;
       end
       function f = get_CP_Function(obj)
           f = @(temp) ((obj.A_Heat_Capacity + obj.B_Heat_Capacity*temp + obj.C_Heat_Capacity*temp^2 + obj.D_Heat_Capacity*temp^3 + obj.E_Heat_Capacity*temp^4)/1000); %(J/(mol*K)
       end
       function setRxn(obj)
           deltaGibbs(obj);
           deltaEnthalpy(obj);
           Delta_Heat_Capacities(obj);
       end
    end
    
    %% Priavte Methods
    methods (Access=private)
        %% Calculate delta G
        function deltaG = deltaGibbs(obj)
            p = size(obj.Product);
            Gibbs_Product = 0;
            for i=1:1:p(2)
                Gibbs_Product = Gibbs_Product + obj.Product(i).Gibbs_Formation;
            end
            r = size(obj.Reactant);
            Gibbs_Reactant = 0;
            for i=1:1:r(2)
                Gibbs_Reactant = Gibbs_Reactant + obj.Reactant(i).Gibbs_Formation;
            end
            obj.Delta_Gibbs = Gibbs_Product - Gibbs_Reactant;
            deltaG = obj.Delta_Gibbs;
        end
        
        %% Calculate delta H
        function deltaH = deltaEnthalpy(obj)
            p = size(obj.Product);
            Enthalpy_Product = 0;
            for i=1:1:p(2)
                Enthalpy_Product = Enthalpy_Product + obj.Product(i).Enthalpy_Formation;
            end
            r = size(obj.Reactant);
            Enthalpy_Reactant = 0;
            for i=1:1:r(2)
                Enthalpy_Reactant = Enthalpy_Reactant + obj.Reactant(i).Enthalpy_Formation;
            end
            obj.Delta_Enthalpy = Enthalpy_Product - Enthalpy_Reactant;
            deltaH = obj.Delta_Enthalpy;
        end
        
        %% Calculate deltaA/B/C/D/E
        function  Delta_Heat_Capacities(obj)
            obj.A_Heat_Capacity = 0;
            obj.B_Heat_Capacity = 0;
            obj.C_Heat_Capacity = 0;
            obj.D_Heat_Capacity = 0;
            obj.E_Heat_Capacity = 0;
            p = size(obj.Product);
            for i=1:1:p(2)
                obj.A_Heat_Capacity = obj.A_Heat_Capacity + obj.Product(i).A_Heat_Capacity;
                obj.B_Heat_Capacity = obj.B_Heat_Capacity + obj.Product(i).B_Heat_Capacity;
                obj.C_Heat_Capacity = obj.C_Heat_Capacity + obj.Product(i).C_Heat_Capacity;
                obj.D_Heat_Capacity = obj.D_Heat_Capacity + obj.Product(i).D_Heat_Capacity;
                obj.E_Heat_Capacity = obj.E_Heat_Capacity + obj.Product(i).E_Heat_Capacity;
                
            end
            r = size(obj.Reactant);
            for i=1:1:r(2)
                obj.A_Heat_Capacity = obj.A_Heat_Capacity - obj.Reactant(i).A_Heat_Capacity;
                obj.B_Heat_Capacity = obj.B_Heat_Capacity - obj.Reactant(i).B_Heat_Capacity;
                obj.C_Heat_Capacity = obj.C_Heat_Capacity - obj.Reactant(i).C_Heat_Capacity;
                obj.D_Heat_Capacity = obj.D_Heat_Capacity - obj.Reactant(i).D_Heat_Capacity;
                obj.E_Heat_Capacity = obj.E_Heat_Capacity - obj.Reactant(i).E_Heat_Capacity;
            end
        end
    end
        
end