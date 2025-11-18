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

f = 0.05; % friction factor

for i = 1:Nxtotal
    Tw = 298; %first guess
    Tw_new = Tw+20;
    index = Nxtotal-i+1;

    %% Coolant Properties 
    kcoolant = 0.134; % coolant thermal conductivity  --> average across experimental data 
    rho_c = 785.26;     % coolant density kg/m^3
    mu_c = 2.038*10^-3; % coolant viscosity 
    %(1/1000)*10^(-0.7009+(841.5/Tc(i))-(Tc(i)*8.6068E-3)-((Tc(i)^2)*8.6068E-3));     % coolant viscosity Pa*s
    %cp_c  = %M*(72.525+ 0.79553*Tc(index+1)-(Tc(index+1)^2)*(2.633E-3)+(3.6498E-6)*Tc(index+1)^3); % coolant heat capacity J/mol K !!!!!!!!!!!!!!!!!!!!!!!
    cp_c = 2000; 
    Pr_c = cp_c * mu_c / kcoolant;
    A_c  = w * d;                      % flow area
    Dh   = 2 * w * d / (w + d);        % hydraulic diameter
    V_c  = mchan / (rho_c * A_c); 



    while abs(Tw_new-Tw) > 15 
         Tw = Tw_new;
         %bartz equation
         sigma(index) = 1/(((0.5*(Tg(index)./T0)*(1+0.5*(k-1)*Mach(index)^2)+0.5)^0.68)*(1+ 0.5*(k-1)*Mach(index)^2)^0.12);
         hg(index) = ((0.026/dt^0.2)*((mu_g(index)^0.2)*cp_g/Pr^0.6)*((Pc*g/Cstar)^0.8)*((dt/R)^0.1))*(((0.5*dt./RADIUS(index)).^2).^0.9)*sigma(index);
         
         Re_c = rho_c * V_c * Dh / mu_c;            %coolant Reynolds number
         %Nu_c = 0.023*(Re_c^0.8)*(Pr_c^0.4);                 %coolant Nusselt number
         % the above correlation may not be accurate for our purposes 
         Nu_c = ((f/8)*(Re_c-1000)*Pr_c)/(1+12.7*((f/8)^0.5)*((Pr_c^2/3)-1));
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
     %Tw = Tw_new;

     %% calculate pressure drop 
end

%% PLOTS 

figure()
plot(xtotal,hg)
ylabel('hg')

figure()
plot(xtotal,u_c)
ylabel('uc')

figure()
hold on
plot(xtotal,Twallinner)
plot(xtotal,Tc(1:i))
plot(xtotal, RADIUS)
hold off
ylabel('Temperature (K)')
legend('Twall','Tcoolant','Nozzle Contour')

figure()
plot(xtotal, Taw)

figure()
plot(xtotal,Q)
ylabel('heat flux')


[MaxWallTemp, idx] = max(Twallinner);
% 
% figure()
% plot(xtotal, u_c)
% hold on
% plot(xtotal, u_c_f)
% legend('uc','ucf')


% Notes:
% 
% * maybe adiabatic wall temp is calculated wrong?
% could be units or area mismatch 
% cool prop 



% T_g*(1+(k-1)/2*M(t)^2) = T_0;