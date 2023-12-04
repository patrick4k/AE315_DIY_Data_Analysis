%% DIY Data Analysis
clear, clc, close all;

%% Load data
U = ["05" "10" "20"];

date_name = [
    "11-29-2023" "parallel"
    "11-30-2023" "perpendicular"
];

complete_data = dictionary();
for i = 1:length(date_name)
    date = date_name(i, 1);
    name = date_name(i, 2);
    data = cell(1,length(U));
    for j = 1:length(U)
        data{j} = table();
        vel = U(j);
        
        % Load tare
        load(sprintf('Data/%s/Sec2_LiftDrag_%s_airfoil_%sms_tare.mat', date, name, vel))
        data{j}.tare = volData;

        % Load actual data
        load(sprintf('Data/%s/Sec2_LiftDrag_%s_airfoil_%sms.mat', date, name, vel))
        data{j}.raw = volData;
    end
    complete_data(name) = {data};
end

clear U data name date date_name vel volData;

%% Analysis
