%% Atmospheric Re-entry Sim
%% Parameters
m = 10000; % mass of object [kg]
v_0 = 50; % initial velocity [m/s]
a_0 = 0; % initial accelertaion [m/s^2]
h_0 = 50000; % initial altitude [m]
d_rocket = 5; % diameter of capsule [m]
A_rocket = pi * (d_rocket/2)^2; % Cross-sectional area [m^2]
rho_0 = 1.225; % Air Density at sea level [kg/m^3]
C_d_rocket = 0.35; % Drag coefficient for rocket
H = 1000000; % Scale Height
G = 6.67e-11; % Gravitational Constant [m^3/kg/s^2]
R_Earth = 6378000; % Radius of Earth [m]
M_Earth = 5.972 * (10^24); % Mass of Earth [kg]
% Heat Transfer
T_amb_0 = 288; % Ambient Temperature at sea level [K]
c_p_shield = 1000; % Heat Capacity of heat shield [J/kg-K] approx.
m_shield = 200; % Mass of heat shield [kg]
h_coefficient = 20; % Heat Transfer Coefficient [W/m^2-K], convetion cooling
epsilon = 0.85; % Emissivity of heat shield surface [dimensionless]
sigma = 5.67e-8; % Stefan-Boltzmann constant [W/m^2-K^4]
eta = 0.5; % Conversion efficiency from mech. energy to heat, approx. 1
A_shield = A_rocket;
% Parachute
C_d_drogue = 0.85; % drag coefficient for drogue
d_drogue = 25; % diameter of drogue [m]
A_drogue = pi * (d_drogue/2)^2; % Area of Drogue chute [m^2]
C_d_parachute = 1.75; % Drag coefficient for parachute
d_parachute = 68; % diameter of parachute [m]
A_parachute = pi * (d_parachute/2)^2; % Area of Parachute [m^2]

%% Time Setup
dt = 0.01; % Time step [s]
t_final = 1300; % Total time [s]
t = 0 : dt : t_final; % Time array

%% Arrays
v = zeros(size(t));
a = zeros(size(t));
h = zeros(size(t));
F_g = zeros(size(t));
F_d = zeros(size(t));
F_net = zeros(size(t));
rho = zeros(size(t));
g_local = zeros(size(t));
Q_in = zeros(size(t)); % heat power input from drag
T_amb = zeros(size(t)); % Ambient temp.
T_shield = zeros(size(t)); % Temp. of heat shield
Q_conv_hist = zeros(size(t));
Q_rad_hist = zeros(size(t));

%% Initial Conditions
h(1) = h_0;
droguechute_deployment = false;
parachute_deployment = false;
t_drogue_start = NaN; % not deployed yet
t_main_start = NaN; % using NaN to represent not deployed yet
t_drogue_duration = 8; % Drogue delpoyment duration [s]
t_deploy_duration = 5; % Parachute deployment duration [s]
terminal_velocity_detected = false;
T_shield(1) = T_amb_0;  % Start at 300 K (room temp)

%% Main Loop
for i = 2 : length(t)
    % Velocity is always negative when falling (down is negative)
    % Drag is always positive when falling (up is positive)

    % Air Density
    rho(i) = rho_0 * exp(-h(i-1) / H);

    %% Drag Force (up is positive)
    if ~droguechute_deployment && h(i-1) < 10000 && v(i-1) < 0
        droguechute_deployment = true;
        t_drogue_start = t(i); % Log time when drogue deploys
    end


    if ~parachute_deployment && droguechute_deployment && h(i-1) < 3000 && v(i-1) < 0
        parachute_deployment = true; % prarchute deploys
        t_main_start = t(i); % Log time when parachute deploys
    end

    % Effective Area & Cd for drag during parachute deployment
    if parachute_deployment == true && droguechute_deployment == true % Both deployed
        t_since_deploy_2 = t(i) - t_main_start; % diff in time since deployment
        if t_since_deploy_2 < t_deploy_duration % time since deployment is less than duration of deployment
            % t_ratio = t_since_deploy / t_deploy_duration; less smooth deployment than using cos             
            t_ratio = 0.5 * (1 - cos(pi * t_since_deploy_2 / t_deploy_duration)); % Smoother transition than linear
            t_ratio = min(t_ratio, 1);

            A_eff = A_drogue + (t_ratio * (A_parachute - A_drogue)); % effective area of parachute over time
            C_d_eff = C_d_drogue + (t_ratio * (C_d_parachute - C_d_drogue)); % effective Cd of parachute over time
        else
            A_eff = A_parachute; % eff area is parachute area once since_deploy = deploy_duration
            C_d_eff = C_d_parachute; % Cd is parachute Cd for same reason
        end
    elseif droguechute_deployment == true && ~parachute_deployment % Drogue only deployed
        t_since_deploy_1 = t(i) - t_drogue_start;
        if t_since_deploy_1 < t_drogue_duration % time since deployment is less than duration of deployment
            t_ratio = 0.5 * (1 - cos(pi * t_since_deploy_1 / t_drogue_duration)); % More smooth transition
            t_ratio = min(t_ratio, 1);

            A_eff = A_rocket + (t_ratio * (A_drogue - A_rocket)); % effective area of parachute over time
            C_d_eff = C_d_rocket + (t_ratio * (C_d_drogue - C_d_rocket)); % effective Cd of parachute over time
        else
            A_eff = A_drogue; % eff area is parachute area once since_deploy = deploy_duration
            C_d_eff = C_d_drogue; % Cd is parachute Cd for same reason
        end

    else % No Parachute deployed
        A_eff = A_rocket;
        C_d_eff = C_d_rocket;
    end
    F_d(i) = 0.5 * C_d_eff * rho(i) * A_eff * v(i-1)^2; % eq. for drag [N]


    % Gravitational Force (down is negative)
    F_g(i) = -(G * M_Earth * m) / (R_Earth + h(i-1))^2;

    % Net Force & Acceleration
    F_net(i) = F_d(i) + F_g(i);
    a(i) = F_net(i) / m;

    % Second diffy eq.
    v(i) = v(i-1) + a(i) * dt; % Change in height (v) [m/s]
    h(i) = h(i-1) + v(i) * dt; % updates height w/ previous height and velocity [m]
   
    % G-Forces
    g_local = (G * M_Earth) ./ (R_Earth + h).^2;
    g_force = abs(a ./ g_local);

    %% Heat Transfer
    % Ambient temperature decreases with altitude approx (simplified):
    T_amb(i) = T_amb_0 * exp(-h(i-1)/7000); % scale height for temperature (7 km typical)
    % Heat gain from drag
    Q_in(i) = eta * F_d(i) * abs(v(i-1));
    % Heat loss from shield
    Q_conv = h_coefficient * A_shield * (T_shield(i-1) - T_amb(i-1)); % Convective cooling
    Q_rad = epsilon * sigma * A_shield * (T_shield(i-1)^4 - T_amb(i-1)^4); % Radiative cooling

    T_shield(i) = T_shield(i-1) + dt/(m_shield*c_p_shield) * (Q_in(i) - Q_conv - Q_rad); % Update shield temp.
    
    %% Terminal Velocity
    v_tol = 0.0003;
    if ~terminal_velocity_detected && i > 5 && abs(v(i) - v(i-1)) < v_tol && abs(v(i-1) - v(i-2)) < v_tol
        % time > 5, difference between prev. and current velocity is less
        % than tolerance, and the step change before, too
        fprintf('Approximate terminal velocity reached at t = %.2f seconds: %.2f m/s\n', t(i), v(i));
        terminal_velocity_detected = true; % makes it so only runs once
    end

    Q_conv_hist(i) = Q_conv;
    Q_rad_hist(i) = Q_rad;

    if h(i) <= 0
        h(i) = 0;
        v(i) = 0;
        a(i) = 0;
        T_shield(i) = T_shield(i-2); % heat doesn't drop to 0
    end

end

disp(parachute_deployment)
disp(t_main_start)
%% Plot
figure
plot(t, h, "b", LineWidth=2)
xlabel('Time (s)')
ylabel('Altitude (m)')
title('Rocket Altitude vs. Time')
xline(t_drogue_start)
xline(t_main_start)
grid on

figure
plot(t, v, "g", LineWidth = 2)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Rocket Velocity vs. Time')
grid on

figure
plot(t, g_force, "magenta")
xlabel('Time (s)')
ylabel('g-force (g)')
title('G-force Experienced During Re-entry')
grid on


figure
plot(t, T_shield - 273.15, "-r", LineWidth=1) % Convert graph to Celsius
xlabel('Time (s)')
ylabel('Temperature (C)')
title('Temperature on re-entry')
grid on

figure;
plot(t, Q_in, 'r', t, Q_conv_hist, 'b', t, Q_rad_hist, 'g');
legend('Heat In', 'Conv Loss', 'Rad Loss');
title('Heat Flow Components');
