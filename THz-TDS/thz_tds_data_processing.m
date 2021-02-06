close all
clear all
%% old datas
% file_names = ["Feb01_006_average.tim","Feb01_007_average.tim","Jan31_009_average.tim","Jan31_007_average.tim"];
%% new datas1
% file_names = ["Jul17_009_average.tim","Jul17_010_average.tim","Jul17_004_average.tim","Jul17_005_average.tim"];
%% NiOx data
% file_names = ["Aug29_002_average.tim","Aug29_003_average.tim","Aug29_005_average.tim","Aug29_006_average.tim"];
%% MAM data
file_names = ["Sep20_001_average.tim","Sep20_002_average.tim","Sep20_003_average.tim"];
%%
files = [];
% file_names: #1 air, #2 glassforair, #3 glassforito, #4 ito
% files(:,8) air, files(:,6) glassforair, files(:,4) glassforito, files(:,2) ito 
% fft_files(:,1) ITO, fft_files(:,2) Glassforito, fft_files(:,3) glassforair, fft_files(:,4) air
sum_phase_air = [];
sum_phase_glass = [];
k = 1;
dg = 1.07*10^-3; % thickness of glass
ds = 130e-7; % thickness of sample in cm
c = 3*10^8;
for i = 1:length(file_names)
    files = cat(2,table2array(importfile(file_names(i))),files); 
end
    N = length(files(:,1));
    time_length = (files(length(files(:,1)),1)-files(1,1))*10^-12;
    freq_step = 1/time_length;
    freq_axis = [];
    freq_axis(1) = 0;
    for n = 2:length(files(:,1))
        freq_axis(n) = freq_axis(n-1) + freq_step;
    end
    while freq_axis(k) < 3*10^12 % 3 THz neden önemli?
        k= k+1;
    end

for j = 1:length(file_names)
    fft_files(:,j) = (fft(files(:,2*j)))/N; % Neden N'e böldük?
end
    freq_axis_plot = freq_axis/10^12; % the unit of x-axis should be THz
abs_fft_files = abs(fft_files(1:k,:)); 
%% Time Domain
    figure
    hold on;
    grid minor;
%     plot(files(:,2));
%     plot(files(:,4));
%     plot(files(:,6));
    plot(files(:,5),files(:,6)); % ??
    plot(files(:,3),files(:,4));
    plot(files(:,1),files(:,2));
%     plot(files(:,5),files(:,8));
    title('Signals in The Time Domain')
    xlabel('Relative Time (ps)')
    ylabel('Electric Field (V/m)')
    %legend('ITO 13\Omega','Glass for ITO','Glass for Air','Air') 
    %legend('NiOx','Glass for NiOx','Glass for Air','Air')
    legend('Air','Quartz','MAM')
%% amplitude
    figure;
    grid minor;
    hold on;
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,3))); 
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,2))); 
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,1))); 
%     plot(freq_axis_plot(1:k),abs(fft_files(1:k,4))); 
    title('Signals in The Frequency Domain');
    ylabel('Amplitude (V/m)')
    xlabel('Frequency (THz)')
    %legend('ITO 13\Omega','Glass for ITO','Glass for Air','Air')
    %legend('NiOx','Glass for NiOx','Glass for Air','Air')  
    legend('Air','Quartz','MAM')
%% Transmittance
    transmittance = ((abs( fft_files(:,1) ) ).^2) ./ (( abs( fft_files(:,2) ) ).^2);
    figure;
    semilogy(freq_axis_plot(1:k),transmittance(1:k));
    grid minor;
    title('Transmittance of The Sample')
    xlim([0 1])
    ylim([1e-3 1e0])
    xlabel('Frequency (THz)')
    ylabel('Transmittance')
%% Absorption Coefficient
    absorp_coeff = (-1/ds).*( log( ((abs( fft_files(:,1) ) ).^2) ./ (( abs( fft_files(:,2) ) ).^2) ) ); 
    % the required thickness should be ds to find the absorption coefficient of the sample?
    figure;
    plot(freq_axis_plot(1:k),absorp_coeff(1:k)); 
    grid minor;
    title('Absorption Coefficient vs Frequency');
    xlim([0 1])
    xlabel('Frequency (THz)')
    ylabel('Absorption Coefficient (1/cm)')
%% Phases for glass
    phase_air = 180*transpose(unwrap(angle(fft_files(:,4))))./pi;
    phase_glass = 180*transpose(unwrap(angle(fft_files(:,3))))./pi;
    figure;
    grid minor;
    hold on;
    plot(freq_axis_plot(1:k),phase_air(1:k));
    hold on;
    plot(freq_axis_plot(1:k),phase_glass(1:k));
    xlim([0 1])
    legend('Phase air','Phase Glass');
    title('Phases');
    ylabel('Phase(\circ)')
    xlabel('Frequency (THz)')
%% Phases for ITO
%     phase_ITO = 180*transpose(phase(fft_files(:,2)))./pi;
%     phase_NiOx = 180*transpose(phase(fft_files(:,1)))./pi;
%     figure;
%     grid minor;
%     hold on;
%     plot(freq_axis_plot(1:k),phase_air(1:k));
%     plot(freq_axis_plot(1:k),phase_glass(1:k));
%     legend('Phase ITO','Phase NiOx');
%     title('phases');
%% Refractive Index of Glass    
    refractive_index = 1+(c*(phase_glass-phase_air)./(2*pi*dg*freq_axis))*pi/180;
    figure;
    plot(freq_axis_plot(1:k),refractive_index(1:k)); 
    xlabel('Frequency (THz)')
    ylabel('Refractive Index')
    xlim([0 1])
    grid minor;
    title('Refractive Index of Glass');
%% Refractive Index of ITO    
%     refractive_index2 = 1+(c*(phase_NiOx-phase_ITO)./(2*pi*(ds+dg)*freq_axis))*3.14/180;
%     figure;
%     plot(freq_axis_plot(1:k),refractive_index2(1:k)); 
%     grid minor;
%     title('refractive index of ITO');
%% Extinction Coefficient
    ext_coeff = ( absorp_coeff ./ (freq_axis.') ) .* (3e10/(4*pi));
    figure
    plot(freq_axis_plot(1:k),ext_coeff(1:k));
    xlim([0 1])
    grid minor
    title('Extinction Coefficient vs. Frequency')
    xlabel('Frequency (THz)')
    ylabel('Extinction Coefficient')
%% sheet conductivity    
    sheet_conductivity = (transpose(refractive_index)+1).*(abs(fft_files(:,2))./abs(fft_files(:,1))-1)/377;
    figure;
    plot(freq_axis_plot(1:k),sheet_conductivity(1:k)); 
    xlim([0 1])
    grid minor;
    xlabel('Frequency (THz)')
    ylabel('Sheet Conductivity (1/(\Omega*sq))')
    title('Sheet Conductivity');    
%% sheet_resistivity
    sheet_resistivity = 1./sheet_conductivity;
    figure;
    plot(freq_axis_plot(1:k),sheet_resistivity(1:k));     
    xlabel('Frequency (THz)')
    xlim([0 1])
    grid minor;
    title('Sheet Resistivity');
    ylabel('Sheet Resistivity (\Omega/sq)')
%% average_sheet_resistivity
    p0=0;
    p1=0;
    for t = 1:length(freq_axis_plot(1,:))
        if gt(freq_axis_plot(t),0.2)
            if p0==0
                p0=t;
            end
        end
        if gt(freq_axis_plot(t),0.8) 
           if p1==0 
                p1=t;
            end 
        end
    end
    sheet_resistivity_sum = 0;
    for c = p0:p1
       sheet_resistivity_sum = sheet_resistivity_sum + sheet_resistivity(c); 
    end
    sheet_resistivity_average = sheet_resistivity_sum/(p1-p0+1);

function average1 = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  AVERAGE1 = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  AVERAGE1 = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  average1 = importfile("C:\Users\atano\Desktop\ITO data\Jan31_20\Jan31_002_average.tim", [7, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 02-Feb-2020 14:08:57

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [7, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["TimeVar", "Var2"];
opts.VariableTypes = ["double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
average1 = readtable(filename, opts);

end