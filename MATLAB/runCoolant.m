function Coolant = runCoolant(Parameter,Geometry,Coolant)
% RUNCOOLANT determine coolant properties for a given temperature and
% pressure.
%
%   Coolant = runCoolant(Parameter,Geometry,Coolant) outputs coolant
%   properties for a given cross section into the Coolant struct.
%
%   See also: RUNHEATTRANSFER

% Index, temperature, and total pressure are stored in local variables for
% convenience.
I = Parameter.index;
T = Coolant.temperature(I); % K
P = Coolant.total_pressure(I); % Pa

%% Density

% Denisty of water and fuel are evaluated using CoolProp.
Water.density(I) = py.CoolProp.CoolProp.PropsSI('D','T', T, 'P', P, 'water');       % kg/m^3
Ethanol.density(I) = py.CoolProp.CoolProp.PropsSI('D','T', T, 'P', P, char(Parameter.fuel));   % kg/m^3
% The denisty of the coolant is determined using a proper formula.
Coolant.density(I) = ((1-Parameter.fuel_ratio)/Water.density(I)+Parameter.fuel_ratio/Ethanol.density(I))^-1; % kg/m^3

% The coolant velocity is determined from the density and geometry.
Coolant.velocity(I) = Parameter.fuel_mass_flow/(Geometry.channel_area(I)*Geometry.channel_number*Coolant.density(I)*cos(Geometry.phi(I))); % m/s

% Static pressure is calculated using Bernoulli's equation.
Coolant.static_pressure(I) = Coolant.total_pressure(I) - 0.5*Coolant.density(I)*Coolant.velocity(I)^2; % Pa

% Static pressure is used for the remaining calculations.
P = Coolant.static_pressure(I); % Pa

%% Water
%  Fluid properties are found for water using CoolProp:
Water.viscosity(I) = py.CoolProp.CoolProp.PropsSI('V','T', T, 'P', P, 'water');                     % Pa-s
Water.conductivity(I) = py.CoolProp.CoolProp.PropsSI('L','T', T, 'P', P, 'water');                  % W/m-K
Water.specific_heat_capacity(I) = py.CoolProp.CoolProp.PropsSI('C','T', T, 'P', P, 'water');        % J/kg-K
Water.surface_tension(I) = py.CoolProp.CoolProp.PropsSI('I','P', P,'Q',0, 'water');                 % N/m
Water.saturation_temperature(I) = py.CoolProp.CoolProp.PropsSI('T', 'P', P, 'Q', 1, 'water');       % K
Water.saturation_pressure(I) = py.CoolProp.CoolProp.PropsSI('P', 'T', T, 'Q', 1, 'water');          % Pa
Water.saturated_liquid_enthalpy(I) = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 0, 'water');    % J/kg
Water.saturated_vapor_enthalpy(I) = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 1, 'water');     % J/kg
Water.saturated_liquid_density(I) = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'Q', 0, 'water');     % kg/m^3
Water.saturated_vapor_density(I) = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'Q', 1, 'water');      % kg/m^3
Water.quality(I) = py.CoolProp.CoolProp.PropsSI('Q','P', P,'T',T, 'water');

%% Ethanol
%  Fluid properties are found for fuel using CoolProp:
Ethanol.viscosity(I) = py.CoolProp.CoolProp.PropsSI('V','T', T, 'P', P, char(Parameter.fuel));                 % Pa-s
Ethanol.conductivity(I) = py.CoolProp.CoolProp.PropsSI('L','T', T, 'P', P, char(Parameter.fuel));              % W/m-K
Ethanol.specific_heat_capacity(I) = py.CoolProp.CoolProp.PropsSI('C','T', T, 'P', P, char(Parameter.fuel));    % J/kg-K
pcrit = py.CoolProp.CoolProp.PropsSI('Pcrit','P', P,'Q',0, char(Parameter.fuel));
if P < pcrit
Ethanol.surface_tension(I) = py.CoolProp.CoolProp.PropsSI('I','P', P,'Q',0, char(Parameter.fuel));             % N/m
Ethanol.saturation_temperature(I) = py.CoolProp.CoolProp.PropsSI('T', 'P', P, 'Q', 1, char(Parameter.fuel));   % K
Ethanol.saturation_pressure(I) = py.CoolProp.CoolProp.PropsSI('P', 'T', T, 'Q', 1, char(Parameter.fuel));      % Pa
Ethanol.saturated_liquid_enthalpy(I) = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 0, char(Parameter.fuel));% J/kg
Ethanol.saturated_vapor_enthalpy(I) = py.CoolProp.CoolProp.PropsSI('H', 'P', P, 'Q', 1, char(Parameter.fuel)); % J/kg
Ethanol.saturated_liquid_density(I) = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'Q', 0, char(Parameter.fuel)); % kg/m^3
Ethanol.saturated_vapor_density(I) = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'Q', 1, char(Parameter.fuel));  % kg/m^3
Ethanol.quality(I) = py.CoolProp.CoolProp.PropsSI('Q','P', P,'T',T, char(Parameter.fuel));
else
Ethanol.surface_tension(I) = py.CoolProp.CoolProp.PropsSI('I','T', T,'Q',0, char(Parameter.fuel));             % N/m
Ethanol.saturation_temperature(I) = py.CoolProp.CoolProp.PropsSI('T', 'P', pcrit, 'Q', 1, char(Parameter.fuel));   % K
Ethanol.saturation_pressure(I) = py.CoolProp.CoolProp.PropsSI('P', 'T', T, 'Q', 1, char(Parameter.fuel));      % Pa
Ethanol.saturated_liquid_enthalpy(I) = py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'Q', 0, char(Parameter.fuel));% J/kg
Ethanol.saturated_vapor_enthalpy(I) = py.CoolProp.CoolProp.PropsSI('H', 'T', T, 'Q', 1, char(Parameter.fuel)); % J/kg
Ethanol.saturated_liquid_density(I) = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', 0, char(Parameter.fuel)); % kg/m^3
Ethanol.saturated_vapor_density(I) = py.CoolProp.CoolProp.PropsSI('D', 'T', T, 'Q', 1, char(Parameter.fuel));  % kg/m^3
Ethanol.quality(I) = py.CoolProp.CoolProp.PropsSI('Q','P', pcrit,'T',T, char(Parameter.fuel));
end

% Throws a warning if the bulk coolant is boiling
if T > Ethanol.saturation_temperature(I)
    warning('Two-phase coolant (boiling)')
end

% Water and fuel fluid properties are combined by mass fraction.
Coolant.viscosity(I) = Ethanol.viscosity(I)*Parameter.fuel_ratio+Water.viscosity(I)*(1-Parameter.fuel_ratio);
Coolant.conductivity(I) = Ethanol.conductivity(I)*Parameter.fuel_ratio+Water.conductivity(I)*(1-Parameter.fuel_ratio);
Coolant.specific_heat_capacity(I) = Ethanol.specific_heat_capacity(I)*Parameter.fuel_ratio+Water.specific_heat_capacity(I)*(1-Parameter.fuel_ratio);
Coolant.surface_tension(I) = Ethanol.surface_tension(I)*Parameter.fuel_ratio+Water.surface_tension(I)*(1-Parameter.fuel_ratio);
Coolant.saturation_temperature(I) = Ethanol.saturation_temperature(I)*Parameter.fuel_ratio+Water.saturation_temperature(I)*(1-Parameter.fuel_ratio);
Coolant.saturation_pressure(I) = Ethanol.saturation_pressure(I)*Parameter.fuel_ratio+Water.saturation_pressure(I)*(1-Parameter.fuel_ratio);
Coolant.saturated_liquid_enthalpy(I) = Ethanol.saturated_liquid_enthalpy(I)*Parameter.fuel_ratio+Water.saturated_liquid_enthalpy(I)*(1-Parameter.fuel_ratio);
Coolant.saturated_vapor_enthalpy(I) = Ethanol.saturated_vapor_enthalpy(I)*Parameter.fuel_ratio+Water.saturated_vapor_enthalpy(I)*(1-Parameter.fuel_ratio);
Coolant.saturated_liquid_density(I) = Ethanol.saturated_liquid_density(I)*Parameter.fuel_ratio+Water.saturated_liquid_density(I)*(1-Parameter.fuel_ratio);
Coolant.saturated_vapor_density(I) = Ethanol.saturated_vapor_density(I)*Parameter.fuel_ratio+Water.saturated_vapor_density(I)*(1-Parameter.fuel_ratio);

end