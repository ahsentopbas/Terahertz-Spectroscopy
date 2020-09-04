close all
clear all
%% old datas
file_names = ["Feb01_010_average.tim","Feb01_009_average.tim","Jan31_008_average.tim","Jan31_007_average.tim"];
%% new datas
% file_names = ["Jul17_009_average.tim","Jul17_010_average.tim","Jul17_004_average.tim","Jul17_005_average.tim"];
%% NiOx data
% file_names = ["Aug29_002_average.tim","Aug29_003_average.tim","Aug29_007_average.tim","Aug29_008_average.tim"];
%%
files = [];
% file_names: #1 air, #2 glassforair, #3 glassforito, #4 ito
% files(:,8) air, files(:,6) glassforair, files(:,4) glassforito, files(:,2) ito 
% fft_files(:,1) ITO, fft_files(:,2) Glassforito, fft_files(:,3) glassforair, fft_files(:,4) air
sum_phase_air = [];
sum_phase_glass = [];
k = 1;
dg = 1.07*10^-3; % thickness of glass
ds = 130e-9; % thickness of sample
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
    plot(files(:,1),files(:,2)); % ??
    plot(files(:,3),files(:,4));
    plot(files(:,5),files(:,6));
    plot(files(:,5),files(:,8));
    legend('ITO','Glassforito','Glassforair','air') 

%% amplitude
    figure;
    grid minor;
    hold on;
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,1))); 
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,2))); 
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,3))); 
    plot(freq_axis_plot(1:k),abs(fft_files(1:k,4))); 
    title('amplitude');
    legend('ITO','GlassforITO','Glassforair','Air') 
%% Absorbance
        absorbance = 1e-3.*(-1/dg).*( log( ((abs( fft_files(:,1) ) ).^2) ./ (( abs( fft_files(:,2) ) ).^2) ) ); 
    figure;
    plot(freq_axis_plot(1:k),absorbance(1:k)); 
    grid minor;
    title('Absorbance vs Frequency');
    xlabel('Frequency (THz)')
    ylabel('Absorbance(?/nm)')
%% Phases for glass
    %phase_air = 180*transpose(phase(fft_files(:,4)))./pi;
    %phase_glass = 180*transpose(phase(fft_files(:,3)))./pi;
        xphase_air = (180/pi).*atan((imag(fft_files(:,4)))./(real(fft_files(:,4))));
        xphase_glass = (180/pi).*atan((imag(fft_files(:,3)))./(real(fft_files(:,3))));
    figure;
    grid minor;
    hold on;
    plot(freq_axis_plot(1:k),xphase_air(1:k));
    plot(freq_axis_plot(1:k),xphase_glass(1:k));
    legend('Phase air','Phase Glass');
    title('Phases');
    ylabel('Phase(\pi)')
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
        refractive_index = 1+(c*(xphase_glass-xphase_air)./(2*pi*dg*freq_axis))*pi/180;
    figure;
    plot(freq_axis_plot(1:k),refractive_index(1:k)); 
    grid minor;
    title('refractive index of glass');
%% Refractive Index of ITO    
%     refractive_index2 = 1+(c*(phase_NiOx-phase_ITO)./(2*pi*(ds+dg)*freq_axis))*3.14/180;
%     figure;
%     plot(freq_axis_plot(1:k),refractive_index2(1:k)); 
%     grid minor;
%     title('refractive index of ITO');
%% Extinction Coefficient
    ext_coeff = ( absorbance ./ (freq_axis.') ) .* (3e8/(4*pi)); %absorbance absorption coefficient de?il, formül yanl??
    figure
    plot(freq_axis_plot(1:k),ext_coeff(1:k));
    grid minor
    title('Extinction Coefficient vs. Frequency')
    xlabel('Frequency (THz)')
    ylabel('Extinction Coefficient')
%% Sheet conductivity    
    sheet_conductivity = (transpose(refractive_index)+1).*(abs(fft_files(:,2))./abs(fft_files(:,1))-1)/377;
    figure;
    plot(freq_axis_plot(1:k),sheet_conductivity(1:k)); 
    grid minor;
    title('sheet conductivity');    
%% sheet_resistivity
    sheet_resistivity = 1./sheet_conductivity;
    figure;
    plot(freq_axis_plot(1:k),sheet_resistivity(1:k));     
    grid minor;
    title('sheet resistivity');
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
    sheet_resistivity_avarage = sheet_resistivity_sum/(p1-p0+1);
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