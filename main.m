%% DIY Data Analysis
clear, clc, close all;

%% Load data

velstr = ["05" "10" "20"];
date_name = [
    "11-29-2023" "parallel"
    "11-30-2023" "perpendicular"
];

% Allocate Space
N = length(date_name) * length(velstr);
U = zeros(N,1);
names = strings(N,1);
Kd_lift = zeros(N,1);
Kd_drag = zeros(N,1);
Vtare_lift = zeros(N,1);
Vtare_drag = zeros(N,1);
Vtare_lift_error = zeros(N,1);
Vtare_drag_error = zeros(N,1);
Vtare_sting_lift = zeros(N,1);
Vtare_sting_drag = zeros(N,1);
Vtare_sting_lift_error = zeros(N,1);
Vtare_sting_drag_error = zeros(N,1);
V_lift = zeros(N,1);
V_drag = zeros(N,1);
V_lift_error = zeros(N,1);
V_drag_error = zeros(N,1);
V_sting_lift = zeros(N,1);
V_sting_drag = zeros(N,1);
V_sting_lift_error = zeros(N,1);
V_sting_drag_error = zeros(N,1);

k = 0;
for i = 1:length(date_name)
    date = date_name(i, 1);
    name = date_name(i, 2);

    % Load Drag Calibration Coeff
    load(sprintf('Data/%s/DragCalibCoeff.mat', date));

    % Load Lift Calibration Coeff
    load(sprintf('Data/%s/LiftCalibCoeff.mat', date));

    for j = 1:length(velstr)
        k = k + 1;
        vel = velstr(j);

        U(k) = str2double(vel);
        names(k) = name;
        Kd_lift(k) = pLift(1);
        Kd_drag(k) = pDrag(1);

        % Load tare
        load(sprintf('Data/%s/Sec2_LiftDrag_%s_airfoil_%sms_tare.mat', date, name, vel))
        avgVolData = mean(volData);
        volData_error = std(volData) / sqrt(30000)*1.96;
        Vtare_lift(k) = avgVolData(1);
        Vtare_drag(k) = avgVolData(2);
        Vtare_lift_error(k) = volData_error(1);
        Vtare_drag_error(k) = volData_error(2);

        % Load tare sting
        load(sprintf('Data/%s/Sec2_LiftDrag_sting_%sms_tare.mat', date, vel))
        avgVolData = mean(volData);
        volData_error = std(volData) / sqrt(30000)*1.96;
        Vtare_sting_lift(k) = avgVolData(1);
        Vtare_sting_drag(k) = avgVolData(2);
        Vtare_sting_lift_error(k) = volData_error(1);
        Vtare_sting_drag_error(k) = volData_error(2);

        % Load actual data
        load(sprintf('Data/%s/Sec2_LiftDrag_%s_airfoil_%sms.mat', date, name, vel))
        avgVolData = mean(volData);
        volData_error = std(volData) / sqrt(30000)*1.96;
        V_lift(k) = avgVolData(1);
        V_drag(k) = avgVolData(2);
        V_lift_error(k) = volData_error(1);
        V_drag_error(k) = volData_error(2);

        % Load actual sting
        load(sprintf('Data/%s/Sec2_LiftDrag_sting_%sms.mat', date, vel))
        avgVolData = mean(volData);
        volData_error = std(volData) / sqrt(30000)*1.96;
        V_sting_lift(k) = avgVolData(1);
        V_sting_drag(k) = avgVolData(2);
        V_sting_lift_error(k) = volData_error(1);
        V_sting_drag_error(k) = volData_error(2);

    end
end

clear data name date date_name vel volData volData_error mass pDrag pLift velstr avgVolData i j k;

%% Analysis

% Constants
mu = 1.716e-5; % dynamic viscosity https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html
Patm = 30.17 * 3386.39; % atmospheric pressure Pa
Tatm = (73-32)*5/9 + 273; % atmospheric temperature K
R = 287; % gas constant
rho = Patm / (R * Tatm); % density kg/m3
c = 175e-3; % chord length m
t = c * 0.12; % thickness m
b = 178.8e-3; % section length m
S = b * c; % characteristic area m2

% Dynamic Pressure
q = rho .* U.^2 ./ 2;

% Drag
D = Kd_drag.*((V_drag - Vtare_drag) - (V_sting_drag - Vtare_sting_drag));

% Drag Coefficient
Cd = D ./ (q .* S);

% Lift
L = Kd_lift.*((V_lift - Vtare_lift) - (V_sting_lift - Vtare_sting_lift));

% Lift Coefficient
Cl = L./(q .* S);

% Reynolds Number
Re = rho .* U .* c ./ mu;


%% Output
fprintf('type\t%s\n', sprintf('%s\t', names));
fprintf('vel\t%s\n', sprintf('%d m/s\t\t', U));
fprintf('D\t%s\n', sprintf('%s\t', D));
fprintf('L\t%s\n', sprintf('%s\t', L));
fprintf('Cd\t%s\n', sprintf('%s\t', Cd));
fprintf('Cl\t%s\n', sprintf('%s\t', Cl));
fprintf('Re\t%s\n', sprintf('%s\t', Re));

%% Plots
close all;
theme = 'dark';

par = find(strcmp(names, "parallel"));
perp = find(strcmp(names, "perpendicular"));

leg = ["Parallel", "Perpendicular"];

% Drag vs U
fig = figure;
fig.Theme = theme;
plot(U(par), D(par), '-.o');
hold on
plot(U(perp), D(perp), '--o');
legend(leg);
grid on;
xlabel('Airspeed [m/s]');
ylabel('Drag [N]');
title('Drag vs Airspeed');

% Drag Ceofficient vs U
fig = figure;
fig.Theme = theme;
plot(U(par), Cd(par), '-.o');
hold on
plot(U(perp), Cd(perp), '--o');
legend(leg);
grid on;
xlabel('Airspeed [m/s]');
ylabel('Cd');
title('Drag Coefficient vs Airspeed');

% Reynolds Number vs U
fig = figure;
fig.Theme = theme;
plot(U(par), Re(par), '-.o');
hold on
plot(U(perp), Re(perp), '--o');
legend(leg);
grid on;
xlabel('Airspeed [m/s]');
ylabel('Re');
title('Reynolds Number vs Airspeed');
