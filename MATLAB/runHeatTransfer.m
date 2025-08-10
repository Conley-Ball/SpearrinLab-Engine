function [Parameter,Geometry,Gas,Coolant,Thermal] = runHeatTransfer(Parameter,Geometry,Gas)
% RUNHEATTRANSFER Simulate the engine during a fire.
%   [Parameter,Geometry,Gas,Coolant,Thermal] = runHeatTransfer(Parameter,Geometry,Gas)
%   simulates the behaviour of the cooling as the engine is being fired.
%   Outputs wall temperatures, coolant properties, and heat/other transport
%   phenomena.
%
%   See also: RUNCOOLANT, BARTZ

% The first index for coolant temperature is taken as ambient.
Coolant.temperature(1) = Parameter.ambient_temperature;

% The first index for coolant pressure is guessed to be 1.2 the chamber 
% pressure.
Coolant.total_pressure(1) = 1.2*Parameter.chamber_pressure;

% The entire simulation runs twice. Once with the guessed initial pressure,
% and once with a corrected pressure.
for i=1:3

% Loops over each station from nozzle exit to injector.
for LCV=1:length(Geometry.axial_distance)
    % Loop control variable is used as an index for current calculation.
    Parameter.index = LCV;

    %% Coolant Properties

    % Function runCoolant is called to get coolant properties at the
    % station's temperature and pressure.
    Coolant = runCoolant(Parameter,Geometry,Coolant);

    %% Critical Heat Flux

    % Temperature difference between stauration and bulk temperature is 
    % found.
    delta_T_sub = (Coolant.saturation_temperature(Parameter.index)-Coolant.temperature(Parameter.index)); % K

    % Temperature difference is converted into English units for critical 
    % heat flux correlation.
    delta_T_sub = delta_T_sub*9/5; % F

    % Coolant velocity is also converted
    velocity = Coolant.velocity(LCV) / (12*0.0254); % ft/s

    % Critical heat flux is calculated using the correlation [2].
    Thermal.critical_heat_flux(LCV) = 0.1003+0.05254*(velocity*delta_T_sub)^0.5; % BTU/in^2-s

    % Critical heat flux is converted back into SI units.
    Thermal.critical_heat_flux(LCV) = Thermal.critical_heat_flux(LCV)*1055.06/0.0254^2; % W/m^2

    % Pressure correction factor is calculated [2].
    F_p = 1.17-8.56e-5*(Coolant.total_pressure(LCV)/6894.76);

    % Pressure correction factor is converted back to SI.
    Thermal.critical_heat_flux(LCV) = Thermal.critical_heat_flux(LCV) * F_p; % W/m^2

    %% Coolant Heat Transfer

    % Reynold's number, Prandtl number for the coolant are calculated by 
    % their definitions.
    Coolant.Reynolds_number(LCV) = Geometry.hydraulic_diameter(LCV)*Coolant.velocity(LCV)*Coolant.density(LCV)/Coolant.viscosity(LCV);
    Coolant.Prandtl_number(LCV) = Coolant.viscosity(LCV)*Coolant.specific_heat_capacity(LCV)/Coolant.conductivity(LCV);

    % Nusselt number is found using the Dittus-Boelter relation [3].
    Coolant.Nusselt_number(LCV) = 0.023*Coolant.Reynolds_number(LCV)^0.8*Coolant.Prandtl_number(LCV)^0.4;

    % Coolant heat transfer coefficient is found from the Nusselt number.
    Coolant.heat_transfer_coefficient_single_phase = Coolant.Nusselt_number(LCV)*Coolant.conductivity(LCV)/Geometry.hydraulic_diameter(LCV);

    % The first term of the nucleate boiling heat transfer coefficient is
    % calculated by the modified Chen correlation [3].
    Coolant.heat_transfer_coefficient_nucleate_boiling = 0.00122*(Coolant.conductivity(LCV)^0.79*Coolant.specific_heat_capacity(LCV)^0.45*Coolant.density(LCV)^0.49*Coolant.surface_tension(LCV)^-0.5*Coolant.viscosity(LCV)^-0.29*(Coolant.saturated_vapor_enthalpy(LCV)-Coolant.saturated_liquid_enthalpy(LCV))^-0.24*Coolant.saturated_vapor_density(LCV)^-0.24);
    
    %% Gas Side Heat Transfer (Bartz relation)

    % Local recovery factor is calculated [4].
    r = Gas.Prandtl_number(LCV)^0.33;

    % Recovery factor is calculated using the local factor [4].
    R = (1+r*(Gas.specific_heat_ratio(LCV)-1)/2*Gas.Mach(LCV)^2)/(1+(Gas.specific_heat_ratio(LCV)-1)/2*Gas.Mach(LCV)^2);

    % Adiabatic wall temperature is found using the total temperature and
    % recover factor [4].
    Gas.adiabatic_wall_temperature(LCV) = Gas.total_temperature(LCV)*R;

    % Initial guess for wall conductivity is made.
    if LCV == 1
        % Wall conductivity is initially guessed to be at ambient conditions.
        [Parameter.k_inner_temperature(1:2),~,~,~,~,~] = materialProperties(Parameter.ambient_temperature,Parameter.material);
        [Parameter.k_inner_channel_temperature(1:2),~,~,~,~,~] = materialProperties(Parameter.ambient_temperature,Parameter.material);
    end

    % Wall conductivity convergence: 
    % Tolerance is set to 1 K and the error is determined by the difference
    % between previous and current outputs.
    tol = 10; % K
    error = tol+1;
    old_temp = 0; % K
    while abs(error) > tol
        % 3 temperatures are calculated by finding the roots of the Bartz 
        % function.
        eq = @(Temp) getfield(Bartz(Parameter,Geometry,Gas,Coolant,Thermal,Temp),'x');
    
        % The solution variables holds the following temperatures: inner wall
        % temperature of the coating near channel and rib, along with wall 
        % temperature next to the channel.
        if LCV == 1
        sln = fsolve(eq,[400,400,400],optimset('Display', 'off', 'TolX',1e-4));
        else
        sln = fsolve(eq,sln,optimset('Display', 'off', 'TolX',1e-4));
        end
        % Thermal struct is updated with temperatures and heats from the Bartz
        % function.
        Thermal = Bartz(Parameter,Geometry,Gas,Coolant,Thermal,sln);
        Thermal.coating_temperature(:,LCV) = sln(1:2)';

        % Error is calculated by the difference between prevous and current
        % temperatures
        error = sln(1) - old_temp;
        old_temp = sln(1);
    
        % Wall conductivity is calculated from the local wall temperature.
        [Parameter.k_inner_temperature(1),~,~,~,~,~] = materialProperties(Thermal.inner_wall_temperature(1,LCV),Parameter.material);
        [Parameter.k_inner_temperature(2),~,~,~,~,~] = materialProperties(Thermal.inner_wall_temperature(2,LCV),Parameter.material);
        [Parameter.k_inner_channel_temperature(1),~,~,~,~,~] = materialProperties(Thermal.inner_temperature(1,LCV),Parameter.material);
        [Parameter.k_inner_channel_temperature(2),~,~,~,~,~] = materialProperties(Thermal.inner_temperature(2,LCV),Parameter.material);
    end

    %% Next Station Properties

    % Darcy friction factor is calculated.
    f = (-1.8*log10((Parameter.roughness/Geometry.hydraulic_diameter(LCV)/3.7)^1.11+6.9/Coolant.Reynolds_number(LCV)))^-2;

    % Expansion or contraction loss coefficient is calculated based on
    % hydraulic diameter
    if not(LCV == length(Geometry.axial_distance))
        % Contraction
        if Geometry.hydraulic_diameter(LCV) > Geometry.hydraulic_diameter(LCV+1)
            beta = Geometry.hydraulic_diameter(LCV)/Geometry.hydraulic_diameter(LCV+1);
            Geometry.loss_coefficient = (beta^2-1)^2;
        % Expansion
        else
            beta = Geometry.hydraulic_diameter(LCV+1)/Geometry.hydraulic_diameter(LCV);
            Geometry.loss_coefficient = 0.5-0.167*(beta)-0.125*(beta)^2-0.208*(beta)^3;
        end
    else
        Geometry.loss_coefficient = 0;
    end

    % Next station coolant pressure is found using friction factor and loss
    % coefficient.
    Coolant.total_pressure(LCV+1) = Coolant.total_pressure(LCV) - 0.5*Coolant.density(LCV)*Coolant.velocity(LCV)^2*(Geometry.loss_coefficient + f*Geometry.dx(LCV)/Geometry.hydraulic_diameter(LCV));

    % Next station coolant temperature is found from the heat flux.
    Coolant.temperature(LCV+1) = Coolant.temperature(LCV) + pi*Geometry.diameter(LCV)*Geometry.dx(LCV)*Thermal.heat_flux(LCV)/(Parameter.fuel_mass_flow*Coolant.specific_heat_capacity(LCV));
end

% Initial coolant pressure is adjusted so that the final pressure matches
% chamber pressure, and the loop is run once more.
Coolant.total_pressure(1) = Coolant.total_pressure(1) + 1.2*Parameter.chamber_pressure - Coolant.total_pressure(end);

% Progress update.
fprintf('.');
end

end

% References:
% [2]  Meyer, Linne, Rouscar. Forced Convection Boiling and Critical Heat Flux of Ethanol in Electrically Heated Tube Tests
% [3] Fang, Yuan, et al. Review of correlations for subcooled flow
% [4] Huzel, Dieter. Huang, David. Modern Engineering for Design of Liquid-Propellant Rocket Engines