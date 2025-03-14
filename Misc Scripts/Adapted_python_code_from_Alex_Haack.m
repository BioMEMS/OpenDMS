% Adapted from Python code by Alexander Haack
% Reference:
% Haack, A., & Hopkins, W. S.: Kinetics in DMS: Modeling Clustering and Declustering Reactions.
% J. Am. Soc. Mass Spectrom. 2022, 33(12), 2250–2262. https://doi.org/10.1021/jasms.2c00224

clc;
clear all;

% Prompt the user to specify the Savitzky-Golay filter level or turn it off
filterOption = input('Enter the Savitzky-Golay filter level (must be an odd number), or enter 0 to skip filtering: ');

% Ask the user if they want to perform background subtraction
bgSubtraction = input('Do you want to perform background subtraction? Enter 1 for Yes, 0 for No: ');

% Specify the folder path containing the .xls files
folderPath = 'C:\Users\Dylan Koch\Box\1_BIOMEMS MEMBERS\KOCH Dylan\0_KOCH\Dylan DMS Paper\Figures\Figure 3';

% Get a list of all .xls files in the folder
fileList = dir(fullfile(folderPath, '*.xls'));

% Loop through each file in the folder
for k = 1:numel(fileList)
    % Create a new figure for each file
    figure('Position', [100, 100, 800, 1000]);
    
    % Construct the full file path
    filePath = fullfile(folderPath, fileList(k).name);
    
    % Get the names of all sheets in the current file
    [~, sheets] = xlsfinfo(filePath);
    
    % Loop through all sheets and process the data
    for i = 1:numel(sheets)
     
        sheetData = readmatrix(filePath, 'Sheet', sheets{i});
        
        yAxisIndices = 1:size(sheetData, 1) - 2; % Indices from 1 to 60

        % Convert indices to the desired range (500 to 1700)
        yAxis = interp1([1, 50], [500, 1700], yAxisIndices);
        y2Axis = interp1([1, 50], [37.18, 148.7], yAxisIndices);

        xData = sheetData(1, 2:end); % X data in the first row
        yDataMatrix = sheetData(3:end, 2:end); % Z data (Ion Current) in the subsequent rows
        
        % Apply background subtraction if selected
        if bgSubtraction
            disp('Background subtraction is being applied.');
            bgSubtractionValue = yDataMatrix(3,10);
            yDataMatrix = yDataMatrix - bgSubtractionValue;
        end
        
        % Apply Savitzky-Golay filtering if selected
        if filterOption > 0 && mod(filterOption, 2) == 1
            yDataMatrix = sgolayfilt(yDataMatrix, 2, filterOption, [], 2);
        elseif filterOption > 0
            disp('Warning: Filter level must be an odd number. Skipping filtering.');
        end
        
        % Create the 3D surface plot using surf
        [X, Y] = meshgrid(xData, yAxis);
        
        surf(X, Y, yDataMatrix);
       
        shading interp; % Optional: smooth shading
        
        colorbar;
        clim([0, max(yDataMatrix(:))]);
        colormap jet;
        
        % Finalize the plot
        title(['Chemical blank: ' fileList(k).name],'FontSize', 20);
        xlabel('Compensation Voltage (CV) [V]','FontSize', 20,'FontWeight','bold');
        ylabel('Separation Voltage (SV) [V]','FontSize', 20,'FontWeight','bold');
        zlabel('Ion Current [pA]','FontSize', 20,'FontWeight','bold');
        ax = gca; % Get current axes
        ax.XAxis.FontSize = 20; % Set the x-axis font size
        ax.YAxis.FontSize = 20; % Set the y-axis font size
        ax.ZAxis.FontSize = 20; % Set the z-axis font size
        view(3); % Set the view to 3D
        xlim([-15 10]);
        ylim([500 1700]);
        view(0,90)
        grid on;

        % Code to compute and plot expected CV vs SVpp values
        % Constants
        kB = 1.380649e-23;  % Boltzmann constant (J/K)

        % Global properties
        T_bath = 319.0;     % Temperature in Kelvin
        p_bath = 101325.0;  % Pressure in Pascals
        Nv = p_bath / (kB * T_bath);  % Particle density (1/m^3)
        gapsize = 0.0005;   % Gap size in meters

        % Molecule mobility values 
         K_lf = 1.933e-4;  % Low-field mobility in m^2/V·s Dmmp
         a_coeffs = [4.0995e-6, -2.3882e-10, +3.3542e-15];  % Alpha coefficients in Td^(-2k)
       

        % Calculate the dispersion plot over a range of SVpp values
        SVpp_grid = linspace(500, 1700, 17);  % SVpp values from 500 to 1700 V in 17 steps
        conv = 0.01;  % Convergence criterion in Volts
        [SVpp_values, CVs] = get_DispPlot(SVpp_grid, conv, K_lf, a_coeffs, Nv, gapsize);

        hold on; % Keep the current figure

        % Set z_value to be maximum Ion Current
        z_value = max(yDataMatrix(:)) * ones(size(CVs));

        % Plot the expected line
        plot3(CVs, SVpp_values, z_value, 'r--', 'LineWidth', 3);

    end
end

%  Function definitions from Alex 
function CV_new = get_CV_it(SVpp, conv, K_lf, a_coeffs, Nv, gapsize)
    % Get CV value (in volts) for a fixed SVpp value (in volts).
    % CVs will be converged to 'conv' volts or until max iterations are reached.
    damp = 0.8;
    maxit = 100;
    % Set up waveform
    Ndt = 1000;
    t_grid = linspace(0, 1, Ndt);
    % Iteration of CVs
    CV_old = 0.0;
    cnt = 0;
    stop = false;
    while ~stop
        E_real = waveform(t_grid, SVpp, gapsize) + CV_old / gapsize;
        K_E = K_func(abs(E_real), K_lf, a_coeffs, Nv);
        vD_av = sum(K_E .* E_real) / Ndt;
        corr = vD_av / K_func(0, K_lf, a_coeffs, Nv) * gapsize;
        CV_new = CV_old - corr * damp;
        cnt = cnt + 1;
        if abs(CV_new - CV_old) < conv || cnt > maxit
            stop = true;
        else
            CV_old = CV_new;
        end
    end
end

function [SVpp_grid, CVs] = get_DispPlot(SVpp_grid, conv, K_lf, a_coeffs, Nv, gapsize)
    Ngrid = length(SVpp_grid);
    CVs = zeros(size(SVpp_grid));
    for i = 1:Ngrid
        CVs(i) = get_CV_it(SVpp_grid(i), conv, K_lf, a_coeffs, Nv, gapsize);
    end
end

function E = waveform(x, SVpp, gapsize)
    % Return field value at certain time, where x=omega*t (dimensionless).
    % x=0-1 is one cycle of the waveform
    D = SVpp * 2 / 3;  % Maximum voltage in sine function in V
    x = x + asin((sqrt(3) - 1) / 2) / (2 * pi);  % Shift by first root, so waveform(0)=0
    f = (2 / 3) * sin(2 * pi * x) + (1 / 3) * sin(4 * pi * x - pi / 2);  % Normalized waveform
    E = D * f / gapsize;  % Electric field in V/m
end

function K_E = K_func(E, K_lf, a_coeffs, Nv)
    % Returns absolute mobility in m^2/V·s at field strength E in V/m
    EN = E ./ Nv * 1e21;  % Reduced field (E/N) in Townsend units (Td)
    K_E = K_lf * ones(size(E));  % Initialize K_E with low-field mobility
    for i = 1:length(a_coeffs)
        K_E = K_E + K_lf .* a_coeffs(i) .* (EN).^(2 * i);
    end
end
