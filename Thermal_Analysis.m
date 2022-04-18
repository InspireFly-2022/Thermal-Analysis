%% InspireFly Thermal Analysis for Low Earth (ISS) Orbit
% Researcher: Nicholas Jones 
% Last Updated: 1/5/2022  
% Topic: Initial Thermal Analysis for a Low Earth Orbit mission
% 
% Resources:  
% STK - Matlab code snippets  
% https://help.agi.com/stkdevkit/index.htm#stkObjects/ObjModMatlabCodeSamples.htm
clear;
clc;
%% Preset Settings  
% Designate if you want the STK Application to open (true) or not (false) 
visibility = input('Enter in if you want the STK application to open (true) or not (false): ');  

% Designate mission length -> true (2 year mission), false (2 day mission)  
fullMission = input('Enter in if you want the mission length to be 2 years (true): '); 
%% Constants 
saveData = true;        % Saves data to .mat and .csv files for later analysis. Set true to enable
makePlots = true;       % Controls whether to make plots of results. Set true to enable
stepSize = 10;          % s - Controls the transient thermal analysis time output step size
exportLength = 1200;    % Controls the number of data points to export for analysis in ANSYS.
                        % This only effects the length of the .csv file, not the .mat file, which will save everything
%% Load In STK Scenario File 
% Create an instance of STK 
stkApplication = actxserver('STK12.Application'); 

% Allows the User to have control over the program 
stkApplication.UserControl = 1; 

% Application Visibility 
if visibility == true  
    stkApplication.Visible = 1; 
else 
    stkApplication.Visible = 0;
end   

% Get the IAgStkObjectRoot interface 
root = stkApplication.Personality2;  

% Load the Scenario   
if fullMission == true 
    root.LoadScenario('C:\Users\William Suffa\Documents\Virginia Tech\Spring 2022\InspireFly\STK Scenarios\InspireFly_Thermal_Analysis.sc'); 
else 
    root.LoadScenario('C:\Users\William Suffa\Documents\Virginia Tech\Spring 2022\InspireFly\STK Scenarios\InspireFly_Thermal_Analysis_2Day.sc') 
end
scenario = root.CurrentScenario;  
 
%1: +X 2: -X 3: +Y 4: -Y 5:+Z 6: -Z
areaF = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01];               % m^2 - Surface area of spacecraft sides (1 -> 6) - 1U CubeSat: 0.01 m^2
epsilonF = [0.96, 0.96, 0.96, 0.96, 0.96, 0.96];                  % Emmisivity of spacecraft surfaces (1 -> 6)
alphaF = [0.96, 0.96, 0.96, 0.96, 0.96, 0.96];                    % Absorptance of spacecraft surfaces (1 -> 6)
surfaceVec = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]; % Normal vectors to spacecraft surfaces defined on body axis.
                                                            % Order: +X, -X, +Y, -Y, +Z, -Z

rho = 1000; % kg m^-3 - spacecraft density
V = 0.001;  % m^3 - spacecraft volume
c = 50;     % J kg^-1 K^-1 - specific heat

% Environment Constants
RE = 6378e3;                % m - Earth equatorial radius
cc = 100;                    % % - Cloud cover percentage of Earth surface
sigma = 5.670374419e-8;     % W m^-2 K^-4 - Stefan-Boltzmann Constant
Tspa = 2.7;                 % K - Temperature of space

% Cloud coverage relationships from:
% Carlson, L., and Horn, W. A New Thermal and Trajectory Model for High Altitude Balloons. 1981.
% Calculate Earth black ball temperature
TE = 214.4 - 0.2 * cc;      % K - Earth black ball temperature

% Calculate Earth albedo
re = 0.18 + 0.0039 * cc;    % albedo coefficent

% Internal heat dissipation from system operation
% Note: currently only support defining as a constant
QOper = 5; % W  
sat = scenario.Children.Item('InspireFly');
root.UnitPreferences.Item('DateFormat').SetCurrentUnit('EpSec');    % Set date format to output in epoch seconds
root.UnitPreferences.Item('Power').SetCurrentUnit('W');             % Set power data to ouput in Watts
    
    % Grab solar intensity data
    sIRptElements = {'Time'; 'Solar Flux'; 'Solar Intensity'};                                      % Define desired report elements. Note that regardless of the order specified here, STK will output the variables alphabetically in the retrieved data set
    satSIDP = sat.DataProviders.GetDataPrvTimeVarFromPath('SEET Vehicle Temperature');              % Get access to data provider in STK
    satSolarData = satSIDP.ExecElements(scenario.StartTime, scenario.StopTime, 1, sIRptElements);   % Have data provider provide desired report elements

    % Grab sun vector data
    sunVecRptElements = {'Time'; 'x/Magnitude'; 'y/Magnitude'; 'z/Magnitude'};
    satSunVecDP = sat.DataProviders.GetDataPrvTimeVarFromPath('Vectors(Body)//Sun');
    satSunVec = satSunVecDP.ExecElements(scenario.StartTime, scenario.StopTime, 1, sunVecRptElements);

    % Grab Earth vector data
    earthVecRptElements = {'Time'; 'x'; 'y'; 'z'};
    satEarthVecDP = sat.DataProviders.GetDataPrvTimeVarFromPath('Vectors(Body)//Nadir(Centric)');
    satEarthVec = satEarthVecDP.ExecElements(scenario.StartTime, scenario.StopTime, 1, earthVecRptElements);

    % Grab subsolar point angle data
    STKEarth = scenario.Children.Item('Earth');   % Create Earth object in scenario
    
    vgtEarth = STKEarth.Vgt;                                % Get reference to vector geometry tool for Earth object
    earth2Sun = vgtEarth.Vectors.Item('Sun(True)');         % Get handle for vecotr from Earth to Sun
    
    vgtSat = sat.Vgt;                                       % Get reference to vector geometry tool for satellite
    sat2Zenith = vgtSat.Vectors.Item('Zenith(Centric)');    % Get handle for vector from satellite to zenith
    
                                                              % Set "To" vector
    
    sSAngleRptElements = {'Time'; 'Angle'};
    satSSAngleDP = sat.DataProviders.GetDataPrvTimeVarFromPath('Angles//AngleFromSubSolarPoint');
    satSSAngle = satSSAngleDP.ExecElements(scenario.StartTime, scenario.StopTime, 1, sSAngleRptElements);
    
    % Grab altitude data
    scAltitudeRptElements = {'Time'; 'Altitude'};
    scAltitudeDP = sat.DataProviders.GetDataPrvTimeVarFromPath('Astrogator Values//Geodetic');
    scAltitude = scAltitudeDP.ExecElements(scenario.StartTime, scenario.StopTime, 1, scAltitudeRptElements);

    % Complete data extraction to usable variables
    STKTime = satSolarData.DataSets.Item(0).GetValues();
    sFlux = satSolarData.DataSets.Item(1).GetValues();
    sIntensity = satSolarData.DataSets.Item(2).GetValues();

    satSunVector = cell(size(STKTime, 1), 3);
    satSunVector( : , 1) = satSunVec.DataSets.Item(1).GetValues();
    satSunVector( : , 2) = satSunVec.DataSets.Item(2).GetValues();
    satSunVector( : , 3) = satSunVec.DataSets.Item(3).GetValues();

    satEarthVector = cell(size(STKTime, 1), 3);
    satEarthVector( : , 1) = satEarthVec.DataSets.Item(1).GetValues();
    satEarthVector( : , 2) = satEarthVec.DataSets.Item(2).GetValues();
    satEarthVector( : , 3) = satEarthVec.DataSets.Item(3).GetValues();

    ssAngle = satSSAngle.DataSets.Item(1).GetValues();
    
    altitude = scAltitude.DataSets.Item(1).GetValues();
% Close STK application upon successful completion of data extraction
clear uiApplication root;

% Clear STK objects
clear sat satEarthVec satEarthVecDP satLightingTimes satSIDP satSSAngle satSSAngleDP
clear satSunVec satSunVecDP scAltitude scAltitudeDP scenario

% Convert from cell arrays to double arrays
STKTime = cell2mat(STKTime);                % s - Time since beginning of scenario
sIntensity = cell2mat(sIntensity);          % % - Percent of Sun disk visible from satellites current position
sFlux = cell2mat(sFlux);                    % W m^-2 - Incoming solar flux
satSunVector = cell2mat(satSunVector);      % Unit vector pointing from satellite to Sun
satEarthVector = cell2mat(satEarthVector);  % Unit vector pointing from satellite to Earth
ssAngle = cell2mat(ssAngle);                % deg - Angle of satellite from subsolar point
hE = cell2mat(altitude) * 10^3;             % m - altitude above Earth surface

n = size(STKTime, 1);                       % Size of STK data arrays

%% Calculate view factors for each surface and Sun angle
% Calculate view factors from spacecraft to Earth
% Reference: http://webserver.dmt.upm.es/~isidoro/tc3/Radiation%20View%20factors.pdf
% "From a small planar surface tilted to a sphere of radius R" scenario

FE = zeros(n, 6);           % View factors for each face to Earth (Surface 1 in column 1 -> Surface 6 in column 6)
sunAngleRaw = zeros(n, 6);  % rad - angle from spacecraft surface normal to Sun vector
sunAngleProc = zeros(n, 6); % rad - angle from spacecraft surface normal to Sun Vector. Angles between 90 and 270 degrees are set to 90 so that they register as zero for the calculations
Beta = zeros(n, 6);         % rad - angle from Spacecraft surface normal vector to Earth vector
y = zeros(n, 6);            % supporting number for view factor calculation - see reference

h = (RE + hE) ./ RE;        % supporting number for view factor calculation - see reference
x = sqrt(h.^2 - 1);         % supporting number for view factor calculation - see reference

for i = 1 : 6
    for j = 1 : n
        Beta( j , i) = acos(dot(surfaceVec(i, : ), satEarthVector(j, :)) / (norm(surfaceVec(i, :)) * norm(satEarthVector(j, :))));
        y ( j , i) = -x(j) * cot(Beta( j, i));
        
        % Test value of beta to see if surface has view to Earth, and then
        % what kind of view
        if abs(pi - Beta(j, i)) <= acos(1 / h(j))
            FE(j, i) = 0;
        elseif abs(Beta(j, i)) <= acos(1 / h(j))
            FE(j, i) = cos(Beta(j, i)) / h(j)^2;
        else
            FE(j, i) = (1 / (pi * h(j)^2)) * (cos(Beta(j, i)) * acos(y(j, i)) - x(j) * sin(Beta(j, i)) * sqrt(1 - y(j, i)^2)) + (1 / pi) * atan(sin(Beta(j, i)) * sqrt(1 - y(j, i)^2) / x(j));
        end
        
        % Also taking advantage of this loop to calculate the angle from
        % each surface normal to the Sun vector
        sunAngleRaw(j, i) = acos(dot(surfaceVec(i, : ), satSunVector(j, : )) / (norm(surfaceVec(i, :)) * norm(satSunVector(j, :))));
        
        % If the angle between the surface normal and the Sun vector is
        % greater than 90 deg (pi / 2 rad) and less than 270 deg (3 pi / 2
        % rad), then the intensity on that surface should be zero. Set sun
        % angle equal to 90 deg (pi / 2 rad) to achieve this
        if sunAngleRaw(j, i) >= pi / 2 && sunAngleRaw(j, i) <= 3 * pi / 2
            sunAngleProc(j, i) = pi / 2;
        else
            sunAngleProc(j, i) = sunAngleRaw(j, i);
        end
    end
    
end

% Calculate view factor from spacecraft to space
% Note: Neglecting view factor to Sun as negligible
Fspa = 1 - FE;

% Plot view factors and angles with time
if makePlots
    % Combined plot for view factors
    figure(1);
    hold on;
    for i = 1 : 6
        plot(STKTime ./ 3600, FE( : , i));
    end
    xlabel('Time (hr)');
    ylabel('F_E');
    title('View Factors to Earth over Analysis Period');
    legend('+X', '-X', '+Y', '-Y', '+Z', '-Z');
    
    % Combined plot for surface beta angles
    figure(2);
    hold on;
    for i = 1 : 6
        plot(STKTime ./ 3600, rad2deg(Beta( : , i)));
    end
    xlabel('Time (hr)');
    ylabel('Surface Beta Angle (deg)');
    title('Surface Beta Angles to Earth over Analysis Period');
    legend('+X', '-X', '+Y', '-Y', '+Z', '-Z');
    
    % Combined plot for angle to Sun
    figure(3);
    hold on;
    for i = 1 : 6
        plot(STKTime ./ 3600, rad2deg(sunAngleRaw( : , i)));
    end
    xlabel('Time (hr)');
    ylabel('Surface Sun Angle (deg)');
    title('Surface Sun Angles over Analysis Period');
    legend('+X', '-X', '+Y', '-Y', '+Z', '-Z');
end

%% Define heat loads
% Expression for Earth IR load
qIR = @(t, Tsc) sigma * (TE^4 - Tsc.^4) .* dot(epsilonF .* areaF, FE(round(t + 1), : ));

% Expression for Earth albedo load
cosdSSAngle = zeros(n, 1);  % Stores cosine of angle from subsolar point for use in albedo function

% Calculate cosine of angle from subsolar point. Places 0 for any angle from
% subsolar point greater than 90 degrees. This would be past the terminator,
% and albedo should no longer be impacting the spacecraft
for i = 1 : n
    if ssAngle(i) > 90
        cosdSSAngle(i) = 0;
    else
        cosdSSAngle(i) = cosd(ssAngle(i));
    end
end

qalbedo = @(t) re * sFlux(round(t + 1)) * (sIntensity(round(t + 1)) ./ 100) .* cosdSSAngle(round(t + 1)) .* dot(alphaF .* areaF, FE(round(t + 1), :));

% Expression for solar load
qSun = @(t) sFlux(round(t + 1)) * (sIntensity(round(t + 1)) ./ 100) .* dot(alphaF .* areaF, cos(sunAngleProc(round(t + 1), :)));

% Expression for heat loss to space
qspa = @(t, Tsc) sigma * (Tspa^4 - Tsc.^4) .* dot(epsilonF .* areaF, Fspa(round(t + 1), : ));

% Build heat balance expression
sigmaq = @(t, Tsc) qIR(t, Tsc) + qalbedo(t) + qSun(t) + qspa(t, Tsc) + QOper;

%% Transient Point Mass Analysis
initialTemp = 295;                                                          % K - Initial temperature of spacecraft
outputTime = STKTime(1) : stepSize : STKTime(end);                          % s - output time vector for transient thermal analysis results
[tResult, TscResult] = ode45(@(t, Tsc) sigmaq(t, Tsc) / (rho * V * c), ...
    outputTime, initialTemp);                                               % Solve ODE for spacecraft temperature using the heat balance equation
m = size(tResult, 1);                                                       % Size of output vectors

%% Extrapolate heat loads versus time
% Find heat loads using result time vector
% Using dot() in the heat loads definitions prevents use of vectorized
% equations
qIRSurf = zeros(m, 6);
qAlbedoSurf = zeros(m, 6);
qSunSurf = zeros(m, 6);
qSpaSurf = zeros(m, 6);

% Calculate loads on each surface as function of time
for k = 1 : m
    for u = 1 : 6
        qIRSurf(k, u) = sigma * (TE^4 - TscResult(k).^4) .* dot(epsilonF(u) .* areaF(u), FE(round(tResult(k) + 1), u));
        qAlbedoSurf(k, u) = re * sFlux(round(tResult(k) + 1)) * (sIntensity(round(tResult(k) + 1)) ./ 100) .* cosdSSAngle(round(tResult(k) + 1)) .* dot(alphaF(u) .* areaF(u), FE(round(tResult(k) + 1), u));
        qSunSurf(k, u) = sFlux(round(tResult(k) + 1)) * (sIntensity(round(tResult(k) + 1)) ./ 100) .* dot(alphaF(u) .* areaF(u), cos(sunAngleProc(round(tResult(k) + 1), u)));
        qSpaSurf(k, u) = sigma * (Tspa^4 - TscResult(k).^4) .* dot(epsilonF(u) .* areaF(u), Fspa(round(tResult(k) + 1), u));
    end
end

% Sum each type of heat load over all surfaces of the spacecraft
qIRSum = sum(qIRSurf, 2);
qAlbedoSum = sum(qAlbedoSurf, 2);
qSunSum = sum(qSunSurf, 2);
qSpaSum = sum(qSpaSurf, 2);

% Sum of heat loads on each surface
surfHeatLoads = qIRSurf + qAlbedoSurf + qSunSurf + qSpaSurf;

%% Plot data
if makePlots
    % Plot temperature
    figure(4);
    plot(tResult ./ (2 * 3600), TscResult);
    grid on;
    xlabel("Time (hr)");
    ylabel("Spacecraft Temperature (K)");
    title("Spacecraft Temperature over 24 Hr Period");
    
    % Combined plot with subplots
    figure(5);
    % Plot temperature
    subplot(2, 1, 1);
    plot(tResult ./ (2*3600), TscResult);
    grid on;
    xlabel("Time (hr)");
    ylabel("Spacecraft Temperature (K)");
    title("Spacecraft Temperature over 24 Hr Period");
    
    % Plot heat loads summed over all surfaces
    subplot(2, 1, 2);
    plot(tResult ./ (2* 3600), qSunSum, tResult ./ (2* 3600), qAlbedoSum, tResult ./ (2* 3600), qIRSum, tResult ./ (2* 3600), qSpaSum);
    grid on;
    xlabel('Time (hr)');
    ylabel('Heat Load (W)');
    title('Heat Load on Spacecraft over 24 Hr Period');
    legend('Solar', 'Albedo', 'Earth IR', 'Space');
    
    % Plot heat loads on each surface
    figure(6);
    hold on;
    for i = 1 : 6
        subplot(6, 1, i);
        plot(tResult ./ (2*3600), qSunSurf( : , i), tResult ./ (2* 3600), qAlbedoSurf( : , i), tResult ./ (2* 3600), qIRSurf( : , i), tResult ./ (2* 3600), qSpaSurf( : , i));
        xlabel('Time (hr)');
        ylabel('Heat Load (W)');
        legend('Solar', 'Albedo', 'IR', 'Space');
        title(strcat('Heat Loads on Surface ', num2str(i)));
    end
    
    % Plot sum of heat loads on each surface
    figure(7);
    hold on;
    for i = 1 : 6
        plot(tResult ./ (2* 3600), surfHeatLoads( : , i));
    end
    xlabel('Time (hr)');
    ylabel('Heat Load (W)');
    title('Heat Load on each Spacecraft Surface over a 24 Hr Period');
    legend('+X', '-X', '+Y', '-Y', '+Z', '-Z');
end

%% Save Data
if saveData
    % Stores resulting temperature, heat loads, and time variable in a .mat file
    % for later analysis
    thermalResults.Tsc = TscResult;
    thermalResults.t = tResult;
    thermalResults.qSunSum = qSunSum;
    thermalResults.qAlbedoSum = qAlbedoSum;
    thermalResults.qIRSum = qIRSum;
    thermalResults.qSpaSum = qSpaSum;
    thermalResults.qSurfaces = surfHeatLoads;
    thermalResults.qSunSurf = qSunSurf;
    thermalResults.qAlbedoSurf = qAlbedoSurf;
    thermalResults.qIRSurf = qIRSurf;
    thermalResults.qSpaSurf = qSpaSurf;
    
    matFileName = strcat(datestr(datetime('now'), 'yyyymmddHHMMSS'), '_Thermal_Analysis_Results.mat');
    save(matFileName, 'thermalResults');
    
    % Save surface heat loads (W) to csv file for analysis in ANSYS
    % tResult has +10 added to it format it for indicating step end times
    exportHeatLoads = [tResult(1 : exportLength) + 10, surfHeatLoads(1 : exportLength, : )];
    csvFileName = strcat(datestr(datetime('now'), 'yyyymmddHHMMSS'), '_Surface_Heat_Loads_100_CC.csv');
    writematrix(exportHeatLoads, csvFileName);
end