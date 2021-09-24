close all;
clc;
clear;

%% Input and result path
Input = '.\Data\Input\';
dir_input = dir([Input,'\*.csv']);
input_number = size(dir_input,1);
Output = '.\Data\Output\';

%% Analysis parameters
% Only take soot mass in the range of (m_min, m_max) for TGA analysis
gaussian_window = 5;    % Parameter to smooth the curve
m_max = 0.9;            % Soot mass limit for analysis
m_min = 0.1;            % Soot mass limit for analysis

%% TGA parameters
heat_rate = 15;  % Heating rate of TGA, K/min
R = 8.314;       % Universal constant, R=8.314J/(mol.K)
n = 1;           % Reaction order of soot oxidation

%% Create a csv file to save the mean data
Oxidation = [Output, '\OxidativeReactivity.csv'];
fid = fopen(Oxidation, 'a+');
% Write header in the csv file
fprintf(fid,'%s,','Case'); % Case name
fprintf(fid,'%s,','Ea');   % Ea
fprintf(fid,'%s,','A');    % A
fprintf(fid,'%s\n','R2');   % R2

for i = 1:input_number
    fn = dir_input(i).name;    % Data name wi. extension
    fn = strrep(fn,'.csv',''); % Data name wo. extension
    ext = '.csv';              % Input data extension
    
    % Results path. If not exist, just create it
    Results = [Output, fn, '\'];
    if exist(Results,'dir') == 0
        mkdir(Results);
    end
    
    % Read the TG data, mass and temperature
    Data_path = [Input, fn, ext];
    TGA_data = textread(Data_path, '', 'delimiter', ',');
    
    Temp = TGA_data(:, 1);   % Temperature in degC
    Temp = Temp + 273.15;    % Convert tempearture from degC into K
    mass = TGA_data(:, 2);   % Unnormalized soot mass
    mass = mass./max(mass);  % Normalize soot mass in range (0, 1)
    
    % Smooth the mass curve
    mass_smooth = smoothdata(mass,'gaussian', gaussian_window); 
    
    % Only soot mass in range of (m_min, m_max) are used for analysis
    x = find(mass_smooth < m_max, 1, 'first');
    y = find(mass_smooth < m_min, 1, 'first');
    mass_smooth = mass_smooth(x:y);
    Temp_smooth = Temp(x:y);    
    
    % Extend the data, the interval of temperature is 0.0001 K
    Temp_ext = Temp_smooth(1):0.0001:Temp_smooth(length(Temp_smooth));
    mass_ext = interp1(Temp_smooth, mass_smooth, Temp_ext, 'spline');
    
    T = Temp_ext';   % Temperature ready for analysis, K
    m = mass_ext';   % Soot mass ready for analysis, (m_min, m_max)
    % Convert temperature into oxidation time in second
    t = (T-T(1)) ./ heat_rate .* 60;  % Oxidation time, s
    
    % Get DTG curve
    dm = diff(m);          % 
    dt = diff(t);          %
    dT = diff(T);
    dm_dt = dm ./ dt;      % dm/dt
    
    m_mid = m(1: length(m)-1 ) + dm;
    t_mid = t(1: length(t)-1 ) + dt;
    T_mid = T(1: length(T)-1 ) + dT;
        
    %% Show the original and smoothed mass curves
    figure()
    
    yyaxis left
    plot(Temp, mass, ':k');
    hold on;
    plot(Temp_smooth, mass_smooth, '--');
    plot(T, m, '-r');
    xlim([min(Temp)-20, max(Temp)+20])
    ylim([-0.05 1.05]);
    xlabel('Temperature (K)');
    ylabel('Soot mass');
    grid on;
    
    yyaxis right
    plot(T_mid, dm_dt, '-g');
    ylim([min(dm_dt)-0.0005 max(dm_dt)+0.0005]);
    ylabel('dm/dt');
    hold off;
    legend('Original', 'Smoothed', 'Extended', 'DTG');
    saveas(gcf, [Results, '1.png']);   
    
    %% Clear useless variables
    clear TGA_data mass mass_smooth Temp Temp_ext mass_ext;   
    
    %% Oxidation rate, alpha
    alpha = (m(1) - m) ./ m(1);

    %% Get X and Y for linear fitting
    Y1 = log(diff(alpha) ./ diff(T));
    Y2 = log((1-alpha) .^ n);
    Y = Y1 - Y2(1:(length(Y2)-1));
    
    X = 1 ./ T;
    X = X(1:(length(X)-1));
    
    %% Linear fitting
    P = polyfit(X, Y, 1);     % Linear regression
    
    k= P(1);
    B = P(2);
    
    % Determination coefficient
    Y_fit = k.* X + B;
    
    SSres = sum( (Y-Y_fit) .^2 );
    SStol = sum( (Y-mean(Y)) .^2 );
    
    R2= 1 - SSres ./ SStol;
   
    % Plot the fitting line
    figure()
    plot(X, Y, '-k');
    hold on;
    plot(X, Y_fit, '-r');
    hold off;
    xlabel('1/T (1/K)');
    ylabel('ln(-dm/dt.m)');
    legend('Original data', 'Linear fit')
    grid on;
    saveas(gcf, [Results, '2.png']);
    
    % Get Ea and A
    Ea = k .* R .* -1e-3;
    A = exp(B) * heat_rate/60;
    
%% Save data
    fprintf(fid, '%s,', fn);
    fprintf(fid, '%4.3f,', Ea);
    fprintf(fid, '%4.3f,', A);
    fprintf(fid, '%4.3f\n', R2);   

    pause(2);
    close all;
    
    % Display how many TG data have been processed.
    X = [num2str(i), '/', num2str(input_number), ' data have been processed.'];
    disp(X);
    
end

fclose all;
disp('Calculation finished.')
