clear all;
close all;
clc;
format long g;

%% Initial data needed

% rover data
lat = 51.0794; lon = -114.133;
roverCoords = [-1641891.178, -3664882.194, 4939966.860];

% rotation matrix
R = [-sind(lon),cosd(lon),0;
    -sind(lat)*cosd(lon), -sind(lat)*sind(lon),cosd(lat);
    cosd(lat)*cosd(lon), cos(lat)*sin(lon), sind(lat)];

% load  output files
errors1 = load('Part1.txt');
errors2 = load('Part2.txt');
errors3 = load('Part3.txt');

%% Task 1 - Time series plot of the single point positioning errors

% calculating  errors in xyz
eX = errors1(:,1) - roverCoords(1,1);
eY = errors1(:,2) - roverCoords(1,2);
eZ = errors1(:,3) - roverCoords(1,3);
errorsXYZ = [eX, eY, eZ];

% calculating errors in enu
errorsENU = (R*(errorsXYZ(:,1:3)'))';
eE = errorsENU(:,1); eN = errorsENU(:,2); eU = errorsENU(:,3);

epochs = 1:size(errors1, 1);

% plotting xyz errors
figure

subplot (3,1,1);
hold on
plot(epochs, eX(:)*10^-3, 'r');
title('Single Point Positioning Error in X direction (km)');
xlabel('Time (s)');
ylabel('Error (km)');

subplot(3,1,2);
plot(epochs, eY(:)*10^-3, 'r');
title('Single Point Positioning Error in Y direction (km)');
xlabel('Time (s)');
ylabel('Error (km)');

subplot(3,1,3);
plot(epochs, eZ(:)*10^-3, 'r');
title('Single Point Positioning Error in Z direction (km)');
xlabel('Time (s)');
ylabel('Error (km)');

hold off

% plotting enu errors
figure

subplot (3,1,1);
hold on
plot(epochs, eE(:)*10^-3, 'r');
title('Single Point Positioning Error in Easting (km) ');
xlabel('Time (s)');
ylabel('Error (km)');

subplot(3,1,2);
plot(epochs, eN(:)*10^-3, 'r');
title('Single Point Positioning Error in Northing (km) ');
xlabel('Time (s)');
ylabel('Error (km)');

subplot(3,1,3);
plot(epochs, eU(:)*10^-3, 'r');
title('Single Point Positioning Error in Vertical (km)');
xlabel('Time (s)');
ylabel('Error (km)');

hold off

%% Task 2

% calculating errors in xyz
eX = errors2(:,1) - roverCoords(1,1);
eY = errors2(:,2) - roverCoords(1,2);
eZ = errors2(:,3) - roverCoords(1,3);
errorsXYZ = [eX, eY, eZ];

% calculatig errors in enu
errorsENU = (R*(errorsXYZ(:,1:3)'))';
eE = errorsENU(:,1);
eN = errorsENU(:,2);
eU = errorsENU(:,3);

% calculating mean errors
mean_EX2 = mean(eX);
mean_EY2 = mean(eY);
mean_EZ2 = mean(eZ);

% calculating root mean squares
rmseX2 = rms(eX);
rmseY2 = rms(eY);
rmseZ2 = rms(eZ);

% plotting xyz errors
figure

subplot (3,1,1);
hold on
plot(epochs, eX(:),'r');
title('Single Difference Positioning Error in X-direction (m)');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(3,1,2);
plot(epochs, eY(:),'r');
title('Single Difference Positioning Error in Y-direction (m)');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(3,1,3);
plot(epochs, eZ(:),'r');
title('Single Difference Positioning Error in Z-direction (m)');
xlabel('Time (s)');
ylabel('Error (m)');

hold off

% plotting enu errors
figure

subplot (3,1,1);
hold on
plot(epochs, eE(:),'r');
title('Single Difference Positioning Error in Easting (m)');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(3,1,2);
plot(epochs, eN(:),'r');
title('Single Difference Positioning Error in Northing (m)');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(3,1,3);
plot(epochs, eU(:),'r');
title('Single Difference Positioning Error in Vertical (m)');
xlabel('Time (s)');
ylabel('Error (m)');

hold off

%% Task 3

% calculating errors in xyz
eX = errors3(:,1) - roverCoords(1,1);
eY = errors3(:,2) - roverCoords(1,2);
eZ = errors3(:,3) - roverCoords(1,3);

% calculating mean errors
mean_EX3 = mean(eX);
mean_EY3 = mean(eY);
mean_EZ3 = mean(eZ);

% calculating root mean squares
rmseX3 = rms(eX);
rmseY3 = rms(eY);
rmseZ3 = rms(eZ);

% plotting xyz errors
figure

subplot (3,1,1);
plot(epochs, eX(:),'r');
title('Single Difference Positioning Error in X-direction (m)');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(3,1,2);
plot(epochs, eY(:),'r');
title('Single Difference Positioning Error in Y-direction (m)');
xlabel('Time (s)');
ylabel('Error (m)');

subplot(3,1,3);
plot(epochs ,eZ(:),'r');
title('Single Difference Positioning Error in Z-direction (m)');
xlabel('Time (s)');
ylabel('Error (m)');
