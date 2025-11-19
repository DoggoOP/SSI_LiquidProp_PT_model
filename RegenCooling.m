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

    %% Coolant Properties - TEMPERATURE DEPENDENT
    Tc_local = Tc(index+1); % Local coolant temperature
    
    % Temperature-dependent properties
    rho_c = 785.26 * (1 - 0.0007*(Tc_local - 298));  
    mu_c = 2.038e-3 * exp(-0.025*(Tc_local - 298));   
    kcoolant = 0.134 * (1 + 0.0002*(Tc_local - 298)); 
    cp_c = 2000 * (1 + 0.0005*(Tc_local - 298));      
    
    Pr_c = cp_c * mu_c / kcoolant;
    A_c  = w * d;                      
    Dh   = 2 * w * d / (w + d);        
    V_c  = mchan / (rho_c * A_c);  
    
    % Store velocity for plotting
    V_array(index) = V_c;  % ADD THIS LINE
    
    % Calculate Reynolds number with updated properties
    Re_c = rho_c * V_c * Dh / mu_c;

    while abs(Tw_new-Tw) > 15 
         Tw = Tw_new;
         %bartz equation
         sigma(index) = 1/(((0.5*(Tg(index)./T0)*(1+0.5*(k-1)*Mach(index)^2)+0.5)^0.68)*(1+ 0.5*(k-1)*Mach(index)^2)^0.12);
         hg(index) = ((0.026/dt^0.2)*((mu_g(index)^0.2)*cp_g/Pr^0.6)*((Pc*g/Cstar)^0.8)*((dt/R)^0.1))*(((0.5*dt./RADIUS(index)).^2).^0.9)*sigma(index);
         
         %Nu_c = 0.023*(Re_c^0.8)*(Pr_c^0.4);                 %coolant Nusselt number
         % the above correlation may not be accurate for our purposes 
         
         % Calculate friction factor based on flow regime
         if Re_c < 2300 % Laminar flow
             f = 64/Re_c;
         else % Turbulent flow - use Colebrook approximation (Haaland equation)
             f = 1/(-1.8*log10((epsilon/Dh/3.7)^1.11 + 6.9/Re_c))^2;
         end
         
         Nu_c = ((f/8)*(Re_c-1000)*Pr_c)/(1+12.7*((f/8)^0.5)*((Pr_c^(2/3))-1));
         u_c(index) = Nu_c*kcoolant/d;                             %coolant convection coefficient
     
         % fin correction
         m = sqrt(2*u_c(index)*w_b(index)/(kwall));                %fin parameter
         e_f = tanh(m*d/w_b(index))/(m*d/w_b(index));      %fin efficiency
         u_c_f(index) = u_c(index)*(w_b(index)+2*e_f*d)/(tw+w_b(index)); %corrected coolant heat trasnfer coefficient

         H(i) = (1/u_c_f(index) + 1/hg(index) + tw/kwall);%corrected bulk convection coefficient
         q_w(index) = (Taw(index)-Tc(index+1))/H(i);
         Tw_new = Taw(index)-q_w(index)/hg(index);
    end
    
     Twallinner(index) = Tw_new;                                 %wall inner temperature
     Q(index) = q_w(index)*pi*d; % heat flow per unit length
     %Tc(i+1) = Taw(i) - q_w(i)*H;
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

fprintf('Average Friction Factor: %.6f\n', mean(f_array));
fprintf('Min Friction Factor: %.6f\n', min(f_array));
fprintf('Max Friction Factor: %.6f\n\n', max(f_array));

fprintf('Total Channel Length: %.6f m\n', xtotal(end));
fprintf('Channel Hydraulic Diameter: %.6f mm\n', Dh*1000);
fprintf('Channel Flow Area: %.6f mm^2\n', A_c*1e6);
fprintf('Average Coolant Velocity: %.6f m/s\n', mean(mchan./(rho_c*A_c)));
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
plot(xtotal,Twallinner)
plot(xtotal,Tc(1:Nxtotal))
plot(xtotal, RADIUS*100) % scale radius for visibility
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
plot(xtotal, f_array)
xlabel('Axial Position (m)')
ylabel('Friction Factor')
title('Friction Factor Distribution')
grid on

[MaxWallTemp, idx] = max(Twallinner);
fprintf('Maximum Wall Temperature: %.6f K at position %.6f m\n', MaxWallTemp, xtotal(idx));