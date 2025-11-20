clear
close all

%% CONSTANTS

T1 = 298; %Coolant Inlet Temperature (K)
De = 44.92E-3; % Nozzle exit diameter (m)
tw = 1E-3; % Chamber thickness (m)
mf = 0.168; % Fuel mass flow rate (kg/s)
R = 4.68E-3; % Radius of curvature at throat (m)
kwall = 140;%Thermal conductivity of wall (W/Km)


T0 = 2530; % Stagnation temp (K) --> taken from RPA
k = 1.24; %RPA value for specific heat ratio
M = 22.6; % molecular weight (g/mol)?? RPA value
Cstar = 1427; %m/s
Pr = 4*k/(9*k-5); % Prandtl number 
Pc = 2.068e+6; % Chamber pressure (Pa)
g = 1; %pretty sure g represents an imperiall unit conversion --> default to 1 here bc we use SI units

d = 1.5E-3; % Channel height (m)
w = 1.5E-3; % Channel width; (m)
dt = 24.36E-3; %Throat diameter (m)
At = 0.25*pi*dt^2; % Throat area (m^2)

Lcyl = 79.97/1000; % Length of cylindrical section of chamber (m)
xchamb = 0:0.001:Lcyl; % indexing across cylindrical section of chamber
Rchamb = 68.14/(2*1000); % radius of chamber (m) 
Achamb = pi*(Rchamb)^2; % Area of chamber (m^2)

Lc = 155.6/1000; % Length from flange to throat (m)
Lconverge = Lc-Lcyl; % Length of converging section of the nozzle (m)
b = 24.89; % degrees --> angle of converging section
Lnozzle = (41.53/1000); % Length of nozzle (m)
 
xnozzle = 0:0.001:Lnozzle;% indexing distance along nozzle axis
rx = NozzleContour(dt,De,17.96,8,length(xnozzle)); % calling NozzleContour function to get the radius as a function of axial length
Aratio = pi*(rx.^2)/At; %area ratio for axial length along nozzle

cp_g = 1.9E3; %heat capacity of gas in chamber J/kgK --> this is an average across the nozzle from RPA, can modify later 

xconverge = 0:0.001:Lconverge; % indexing distance along converging section
Rconverging = Rchamb - ((Rchamb-dt/2)/Lconverge)*(xconverge); % Radius of converging section as a function of axial length
Aratioconv = (pi*Rconverging.^2)/At; % Area of converging section as a function of axial length

figure()
plot(xnozzle,rx) % this just plots a contour of the nozzle from throat to end

Arearelation = [];
Mach = [];

% this for loop calculates mach number and area in cylindrical section of
% combustion chamber --> we assume M=0 in chamber and Area is constant
for i = 1:length(xchamb)
    if xchamb(i) <= Lcyl 
        Mach(i)=0;
        Arearelation(i)=Achamb/At;
        RADIUS(i)= Rchamb;
    end
end

% this for loop calculates mach number and area in converging (subsonic) section of
% nozzle --> we assume isentropic flow 
for i = 1:length(xconverge)
    if xconverge(i)+Lcyl <= Lc
        [mach, T, P, rho, area] = flowisentropic(k, Aratioconv(i),'sub');
        Arearelation(i+length(xchamb)) = area;
        Mach(i+length(xchamb)) = mach;
        RADIUS(i+length(xchamb)) = Rconverging(i);
    end
end

% this for loop calculates mach number and area in diverging (supersonic) section of
% nozzle --> we assume isentropic flow 
for i = 1:length(xnozzle)
    if xnozzle(i)+Lcyl+Lconverge <= Lconverge+Lcyl+Lnozzle
        [mach, T, P, rho, area] = flowisentropic(k, Aratio(i),'sup');
        Arearelation(i+length(xchamb)+length(xconverge)) = area;
        Mach(i+length(xchamb)+length(xconverge)) = mach;
        RADIUS(i+length(xchamb)+length(xconverge)) = rx(i);
    end
end

A=Arearelation/At;

%% Plot of entire chamber contour 
xtotal = 0:0.001:Lnozzle+Lcyl+Lconverge;
figure()
plot(xchamb,RADIUS(1:length(xchamb)))
hold on
plot(xconverge+xchamb(end),RADIUS(length(xchamb):length(xchamb)+length(xconverge)-1))
hold on
plot(xnozzle+xconverge(end)+xchamb(end),RADIUS(length(xchamb)+length(xconverge):length(xchamb)+length(xconverge)+length(xnozzle)-1))
hold on
plot(xtotal,xtotal*0,'LineStyle','--')
axis equal

%% REGEN ITERATIVE SECTION

% T0 = T0 *(1+ (((k-1)/2))*Mach.^2).^(-1); %mixture temp
Tg = T0./(1+((k-1)/2)*Mach.^2); %Static Temperature of gas as a function of mach number (K)
mu_g = (46.6E-10)*sqrt(M).*Tg.^0.6; %viscosity
%mu_g = 1.716e-5 * (Tg/273.15).^0.7; % viscosity of gas as a function of temperature huzel & huang 
Taw = T0.*((1+Pr^(1/3)*(k-1)./2*Mach.^2)./(1+(k-1)./2*Mach.^2)); %adiabatic wall temp (K)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Tc = T1+55/2; % coolant bulk temp taken as average between inlet and outlet conditions 
% % --> Huzel & Huang gave typically ~55º temp difference 

N = 35; % number of channels 
%N = pi*(De+0.8*(d+2*tw))/(d+2*tw); %  of channels 
mchan = mf/N; % mass flow rate per channel (kg/s)
Perimeter_c = 4*d;         % wetted perimeter per channel (M)
dx = xnozzle(2) - xnozzle(1);          % axial step along nozzle dx

% Tc = 298; %coolant inlet temp
% Tc(1) = Tc;
Nxtotal = numel(xtotal); %number of axial steps
Tc = 298*ones(1,Nxtotal+1); %coolant temp starts at 298K 

% this loop calculates the fin base widht as a function of distance along
% nozzle axis 
for i = 1:length(RADIUS)
    d_i = 2*RADIUS(i);                       %inner diameter
    d_o = d_i + 2*tw;                 %outer diameter
    w_b(i) = (pi*d_o - N*w)/N;   %fin base width
end

% Surface roughness for copper (from hand calcs)
epsilon = 0.002E-3; % 0.002 mm in meters

% NEW APPROACH: Set the REQUIRED injector pressure
% This is the pressure needed at the chamber inlet where coolant exits to the injector
P_injector_required = 0.15e6; % 0.15 MPa = 150 kPa = 21.76 psi (example value, adjust as needed)
Pc_coolant = zeros(1, Nxtotal+1);

for i = 1:Nxtotal
    Tw = 298; %first guess
    Tw_new = Tw+20;
    index = Nxtotal-i+1;

    %% Coolant Properties - TEMPERATURE DEPENDENT FOR IPA
    Tc_local = Tc(index+1); % Local coolant temperature
    
    % IPA Temperature-dependent properties
    rho_c = 1074.4 - 0.9226*Tc_local - 0.00104*Tc_local^2;  
    % Vogel-Fulcher-Tammann equation for IPA viscosity
    A = -3.476;
    B = 726.8;
    C = 133.2;
    mu_c = exp(A + B/(Tc_local - C)) / 1000; % Convert to Pa·s  
    kcoolant = 0.169 - 0.000118*Tc_local; 
    cp_c = 1873 + 4.04*Tc_local;      
    
    Pr_c = cp_c * mu_c / kcoolant;
    A_c  = w * d;                      
    Dh   = 2 * w * d / (w + d);        
    V_c  = mchan / (rho_c * A_c);  
    
    % Store velocity for plotting
    V_array(index) = V_c;
    
    % Calculate Reynolds number with updated properties
    Re_c = rho_c * V_c * Dh / mu_c;

    while abs(Tw_new-Tw) > 15 
         Tw = Tw_new;
         %bartz equation
         sigma(index) = 1/(((0.5*(Tg(index)./T0)*(1+0.5*(k-1)*Mach(index)^2)+0.5)^0.68)*(1+ 0.5*(k-1)*Mach(index)^2)^0.12);
         hg(index) = ((0.026/dt^0.2)*((mu_g(index)^0.2)*cp_g/Pr^0.6)*((Pc*g/Cstar)^0.8)*((dt/R)^0.1))*(((0.5*dt./RADIUS(index)).^2).^0.9)*sigma(index);
         
         % Calculate friction factor based on flow regime
         if Re_c < 2300 % Laminar flow
             f = 64/Re_c;
         else % Turbulent flow - use Colebrook approximation (Haaland equation)
             f = 1/(-1.8*log10((epsilon/Dh/3.7)^1.11 + 6.9/Re_c))^2;
         end
         
         % Calculate Nusselt number based on flow regime - CORRECTED
         if Re_c < 2300 % Laminar flow
             % For rectangular channels, use Shah correlation
             % For developing laminar flow with constant heat flux
             alpha = min(d,w)/max(d,w); % aspect ratio (always <= 1)
             
             % Shah correlation for fully developed laminar flow in rectangular ducts
             if alpha == 1 % Square channel
                 Nu_laminar_fd = 3.608; % Fully developed Nu for square duct
             else
                 % Interpolation for rectangular ducts
                 Nu_laminar_fd = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - ...
                                5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
             end
             
             % Account for thermal entry length effects
             % For developing flow, Nu can be higher
             % Gz = (D_h/x) * Re * Pr (Graetz number)
             x_thermal = i * dx; % distance from entrance
             if x_thermal < 0.001
                 x_thermal = 0.001; % Avoid division by zero
             end
             Gz = (Dh / x_thermal) * Re_c * Pr_c;
             
             if Gz > 100 % Developing flow
                 % Use Hausen correlation for entry length
                 Nu_c = 3.66 + (0.0668 * Gz) / (1 + 0.04 * Gz^(2/3));
             else % Fully developed
                 Nu_c = Nu_laminar_fd;
             end
             
             % Ensure minimum Nusselt number
             if Nu_c < 3.0
                 Nu_c = 3.0;
             end
             
         else % Turbulent flow (Re >= 2300)
             % Gnielinski correlation for turbulent flow
             Nu_c = ((f/8)*(Re_c-1000)*Pr_c)/(1+12.7*sqrt(f/8)*((Pr_c^(2/3))-1));
             
             % For transition region (2300 < Re < 4000), blend between laminar and turbulent
             if Re_c < 4000
                 % Calculate laminar Nu at Re=2300
                 alpha = min(d,w)/max(d,w);
                 if alpha == 1
                     Nu_laminar = 3.608;
                 else
                     Nu_laminar = 7.541 * (1 - 2.610*alpha + 4.970*alpha^2 - ...
                                    5.119*alpha^3 + 2.702*alpha^4 - 0.548*alpha^5);
                 end
                 
                 % Linear interpolation in transition region
                 transition_factor = (Re_c - 2300) / (4000 - 2300);
                 Nu_c = Nu_laminar * (1 - transition_factor) + Nu_c * transition_factor;
             end
         end
         
         u_c(index) = Nu_c*kcoolant/Dh;  % Use Dh instead of d for consistency
     
         % fin correction
         m = sqrt(2*u_c(index)*w_b(index)/(kwall));                %fin parameter
         e_f = tanh(m*d/w_b(index))/(m*d/w_b(index));      %fin efficiency
         u_c_f(index) = u_c(index)*(w_b(index)+2*e_f*d)/(tw+w_b(index)); %corrected coolant heat trasnfer coefficient

         H(i) = (1/u_c_f(index) + 1/hg(index) + tw/kwall);%corrected bulk convection coefficient
         q_w(index) = (Taw(index)-Tc(index+1))/H(i);
         Tw_new = Taw(index)-q_w(index)/hg(index);
         
         % Store Nu for diagnostics
         Nu_array(index) = Nu_c;
    end
    
     Twallinner(index) = Tw_new;                                 %wall inner temperature
     Q(index) = q_w(index)*pi*d; % heat flow per unit length
     qprime = q_w(index) * N * Perimeter_c;        % W/m (total heat to coolant per unit length)
     dT_c(index)   = qprime * dx / (mf * cp_c);       % mf is total coolant mass flow
     Tc(index) = Tc(index+1) + dT_c(index);
     
     %% Calculate pressure drop using friction factor from above
     % Darcy-Weisbach equation
     dP_dx = f * (1/Dh) * (rho_c * V_c^2 / 2);  % Pa/m
     dP(index) = dP_dx * dx;                     % Pressure drop over this segment (Pa)
     f_array(index) = f;                         % Store friction factor for plotting
     Re_array(index) = Re_c;                     % Store Reynolds number
     
     % Don't accumulate pressure yet - we'll do that after the loop
end

%% Now accumulate pressure starting from the required injector pressure
% Chamber inlet (position 1) is where coolant exits to injector
Pc_coolant(1) = P_injector_required;

% Work backwards along flow direction to find pressure at each point
% Pressure increases going upstream (against flow direction)
for index = 2:Nxtotal+1
    Pc_coolant(index) = Pc_coolant(index-1) + dP(index-1);
end

% Total pressure drop
Delta_P_total = Pc_coolant(Nxtotal+1) - Pc_coolant(1);
fprintf('\n========== PRESSURE DROP RESULTS ==========\n');
fprintf('Total Pressure Drop:\n');
fprintf('  %.6f MPa\n', Delta_P_total/1e6);
fprintf('  %.6f Pa\n', Delta_P_total);
fprintf('  %.6f psi\n', Delta_P_total/6894.76);
fprintf('  %.6f bar\n\n', Delta_P_total/1e5);

fprintf('REQUIRED Injector Pressure (Chamber Inlet, Coolant Exit):\n');
fprintf('  %.6f MPa\n', Pc_coolant(1)/1e6);
fprintf('  %.6f Pa\n', Pc_coolant(1));
fprintf('  %.6f psi\n', Pc_coolant(1)/6894.76);
fprintf('  %.6f bar\n\n', Pc_coolant(1)/1e5);

fprintf('REQUIRED Feed Pressure (Nozzle Exit, Coolant Inlet):\n');
fprintf('  %.6f MPa\n', Pc_coolant(Nxtotal+1)/1e6);
fprintf('  %.6f Pa\n', Pc_coolant(Nxtotal+1));
fprintf('  %.6f psi\n', Pc_coolant(Nxtotal+1)/6894.76);
fprintf('  %.6f bar\n\n', Pc_coolant(Nxtotal+1)/1e5);

fprintf('Average Reynolds Number: %.6f\n', mean(Re_array));
fprintf('Min Reynolds Number: %.6f\n', min(Re_array));
fprintf('Max Reynolds Number: %.6f\n\n', max(Re_array));

fprintf('Average Nusselt Number: %.6f\n', mean(Nu_array));
fprintf('Min Nusselt Number: %.6f\n', min(Nu_array));
fprintf('Max Nusselt Number: %.6f\n\n', max(Nu_array));

fprintf('Average Friction Factor: %.6f\n', mean(f_array));
fprintf('Min Friction Factor: %.6f\n', min(f_array));
fprintf('Max Friction Factor: %.6f\n\n', max(f_array));

fprintf('Total Channel Length: %.6f m\n', xtotal(end));
fprintf('Channel Hydraulic Diameter: %.6f mm\n', Dh*1000);
fprintf('Channel Flow Area: %.6f mm^2\n', A_c*1e6);
fprintf('==========================================\n\n');

fprintf('DESIGN SUMMARY:\n');
fprintf('If your injector requires %.2f psi,\n', Pc_coolant(1)/6894.76);
fprintf('your pump/tank must provide at least %.2f psi\n', Pc_coolant(Nxtotal+1)/6894.76);
fprintf('to overcome the %.2f psi friction loss in the cooling channels.\n', Delta_P_total/6894.76);
fprintf('==========================================\n');

%% PLOTS 

figure()
plot(xtotal,hg)
ylabel('hg (W/m^2K)')
xlabel('Axial Position (m)')
title('Gas Side Heat Transfer Coefficient')
grid on

figure()
plot(xtotal,u_c)
ylabel('u_c (W/m^2K)')
xlabel('Axial Position (m)')
title('Coolant Heat Transfer Coefficient')
grid on

figure()
hold on
plot(xtotal,Twallinner, 'LineWidth', 2)
plot(xtotal,Tc(1:Nxtotal), 'LineWidth', 2)
plot(xtotal, RADIUS*100, 'LineWidth', 1) % scale radius for visibility
hold off
ylabel('Temperature (K)')
xlabel('Axial Position (m)')
legend('Twall','Tcoolant','Nozzle Contour (scaled x100)')
title('Temperature Distribution')
grid on

figure()
plot(xtotal, Taw)
ylabel('Adiabatic Wall Temperature (K)')
xlabel('Axial Position (m)')
title('Adiabatic Wall Temperature')
grid on

figure()
plot(xtotal,Q)
ylabel('Heat Flux (W/m)')
xlabel('Axial Position (m)')
title('Heat Flux Distribution')
grid on

% Pressure drop plots
figure()
plot(xtotal, Pc_coolant(1:end-1)/1e6)
xlabel('Axial Position (m)')
ylabel('Coolant Pressure (MPa)')
title('Coolant Channel Pressure Distribution')
grid on

figure()
plot(xtotal, dP/1000)
xlabel('Axial Position (m)')
ylabel('Pressure Drop per Segment (kPa)')
title('Local Pressure Drop Distribution')
grid on

figure()
yyaxis left
plot(xtotal, Pc_coolant(1:end-1)/1e6, 'LineWidth', 2)
ylabel('Coolant Pressure (MPa)', 'FontSize', 12)
xlabel('Axial Position (m)', 'FontSize', 12)
ylim([min(Pc_coolant(1:end-1)/1e6)*0.95, max(Pc_coolant(1:end-1)/1e6)*1.05])
ax = gca;
ax.YColor = 'b';

yyaxis right
plot(xtotal, V_array, 'LineWidth', 2, 'Color', 'r')
ylabel('Coolant Velocity (m/s)', 'FontSize', 12)
ylim([min(V_array)*0.95, max(V_array)*1.05])
ax = gca;
ax.YColor = 'r';

title('Coolant Pressure and Velocity Distribution', 'FontSize', 14)
legend('Pressure', 'Velocity', 'Location', 'best')
grid on

% Additional diagnostic plots
figure()
plot(xtotal, Re_array)
xlabel('Axial Position (m)')
ylabel('Reynolds Number')
title('Reynolds Number Distribution')
grid on

figure()
plot(xtotal, Nu_array, 'LineWidth', 2)
xlabel('Axial Position (m)')
ylabel('Nusselt Number')
title('Nusselt Number Distribution')
grid on
hold on
yline(3.608, 'r--', 'Laminar Fully Developed', 'LineWidth', 1.5)
hold off

figure()
plot(xtotal, f_array)
xlabel('Axial Position (m)')
ylabel('Friction Factor')
title('Friction Factor Distribution')
grid on

[MaxWallTemp, idx] = max(Twallinner);

% DIAGNOSTIC PLOT
figure()
subplot(4,1,1)
plot(xtotal, Taw, 'LineWidth', 2)
ylabel('T_{aw} (K)')
title('Adiabatic Wall Temperature')
grid on

subplot(4,1,2)
plot(xtotal, hg, 'LineWidth', 2)
ylabel('h_g (W/m^2K)')
title('Gas Side Heat Transfer Coefficient')
grid on

subplot(4,1,3)
plot(xtotal, Tc(1:Nxtotal), 'LineWidth', 2)
ylabel('T_c (K)')
title('Coolant Bulk Temperature')
grid on

subplot(4,1,4)
plot(xtotal, Twallinner, 'LineWidth', 2)
hold on
plot(xtotal, Taw, '--')
ylabel('Temperature (K)')
xlabel('Axial Position (m)')
title('Wall and Adiabatic Temperatures')
legend('T_{wall}', 'T_{aw}')
grid on

% Add vertical line at throat
throat_position = Lcyl + Lconverge;
for i = 1:4
    subplot(4,1,i)
    hold on
    xline(throat_position, 'r--', 'Throat', 'LineWidth', 1.5)
end

fprintf('Maximum Wall Temperature: %.6f K at position %.6f m\n', MaxWallTemp, xtotal(idx));