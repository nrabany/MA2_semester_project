% Script written to treat and analyze the data that were recovered with an
% IMU in the context of a semester project made at the EPFL
% Script created by Nicolas Rabany, EPFL student
% Last modified the 08/06/2017
%% Load data
clear all
clc

plot_index = 1; % 0 -> no plots 
                % 1 -> general plots (general, filtered, angles, without
                %       gravity)
                % 2 -> all plots (including details plots for filtering)
% LOAD DATA
A = load('walking.txt'); % put the file that you want analyze
time = A(:,1)/1000;
a_x = -A(:,2); %  /!\ imu upside down
a_y = A(:,3);
a_z = -A(:,4); %  /!\ imu upside down 
omega_x =  -A(:,5); %  /!\ imu upside down
omega_y =  A(:,6);
omega_z = A(:,7); % /!\ imu upside down

if(plot_index == 1 || plot_index ==2)
    % Plot data for general overview
    figure
    % Plot acceleromter data
    subplot(2,1,1);
    hold on
    plot(time, a_x, 'r');
    plot(time, a_y, 'g');
    plot(time, a_z, 'b');
    title('Accelerometer'); xlabel('time(s)'); ylabel('Acceleration(mg)');
    legend('x','y','z');
    hold off
    % Plot gyroscope data
    subplot(2,1,2);
    hold on
    plot(time, omega_x, 'r');
    plot(time, omega_y, 'g');
    plot(time, omega_z, 'b');
    title('Gyroscope'); xlabel('Time(s)'); ylabel('Angular rate(°/s)');
    legend('x','y','z');
end
%%
% FILTER DATA
% Accerometer data filtering
a = 1;
coeff_acc = [0.2 0.6 0.2]; % choose odd number
window_size = length(coeff_acc);
loop = 2;

%filter data with a window of window_size
filt_a_x = filter(coeff_acc, a, a_x);
filt_a_y = filter(coeff_acc, a, a_y);
filt_a_z = filter(coeff_acc, a, a_z);

filt_a_x(1:(end-(window_size-1)/2)) = filt_a_x((window_size-1)/2+1:end);
filt_a_x((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_a_x(end);
filt_a_y(1:(end-(window_size-1)/2)) = filt_a_y((window_size-1)/2+1:end);
filt_a_y((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_a_y(end);
filt_a_z(1:(end-(window_size-1)/2)) = filt_a_z((window_size-1)/2+1:end);
filt_a_z((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_a_z(end);

filt_a_x2 = filt_a_x;
filt_a_y2 = filt_a_y;
filt_a_z2 = filt_a_z;

for i = 1:loop
    filt_a_x2 = filter(coeff_acc, a, filt_a_x2);
    filt_a_y2 = filter(coeff_acc, a, filt_a_y2);
    filt_a_z2 = filter(coeff_acc, a, filt_a_z2);

    filt_a_x2(1:(end-(window_size-1)/2)) = filt_a_x2((window_size-1)/2+1:end);
    filt_a_x2((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_a_x2(end);
    filt_a_y2(1:(end-(window_size-1)/2)) = filt_a_y2((window_size-1)/2+1:end);
    filt_a_y2((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_a_y2(end);
    filt_a_z2(1:(end-(window_size-1)/2)) = filt_a_z2((window_size-1)/2+1:end);
    filt_a_z2((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_a_z2(end);

end

% Gyroscope data filtering
a2 = 1;
coeff_gyr = [0.2 0.6 0.2]; % chose odd number
window_size2 = length(coeff_gyr);

%filter data with a window of window_size
filt_omega_x = filter(coeff_gyr, a2, omega_x);
filt_omega_y = filter(coeff_gyr, a2, omega_y);
filt_omega_z = filter(coeff_gyr, a2, omega_z);

filt_omega_x(1:(end-(window_size-1)/2)) = filt_omega_x((window_size-1)/2+1:end);
filt_omega_x((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_omega_x(end);
filt_omega_y(1:(end-(window_size-1)/2)) = filt_omega_y((window_size-1)/2+1:end);
filt_omega_y((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_omega_y(end);
filt_omega_z(1:(end-(window_size-1)/2)) = filt_omega_z((window_size-1)/2+1:end);
filt_omega_z((end-(window_size-1)/2+1):end) = ones(1,(window_size-1)/2)*filt_omega_z(end);

if(plot_index == 1 || plot_index ==2)
    % Plot filtered data
    figure
    % Plot acceleromter data
    subplot(2,1,1);
    hold on
    plot(time, filt_a_x2, 'r');
    plot(time, filt_a_y2, 'g');
    plot(time, filt_a_z2, 'b');
    title('Filtered Accelerometer'); xlabel('time(s)'); ylabel('Acceleration(mg)');
    legend('x','y','z');
    hold off
    % Plot gyroscope data
    subplot(2,1,2);
    hold on
    plot(time, filt_omega_x, 'r');
    plot(time, filt_omega_y, 'g');
    plot(time, filt_omega_z, 'b');
    title('Filtered Gyroscope'); xlabel('Time(s)'); ylabel('Angular rate(°/s)');
    legend('x','y','z');
    hold off
end

if(plot_index == 2)
    figure
    subplot(3,2,1);
    plot(time,a_x)
    hold on
    plot(time,filt_a_x2)
    legend('a_x','filtered a_x');
    subplot(3,2,3);
    plot(time,a_y)
    hold on
    plot(time,filt_a_y2)
    legend('a_y','filtered a_y');
    subplot(3,2,5);
    plot(time,a_z)
    hold on
    plot(time,filt_a_z2)
    legend('a_z','filtered a_z');
    subplot(3,2,2);
    plot(time,omega_x)
    hold on
    plot(time,filt_omega_x)
    legend('\omega_x','filtered \omega_x');
    subplot(3,2,4);
    plot(time,omega_y)
    hold on
    plot(time,filt_omega_y)
    legend('\omega_y','filtered \omega_y');
    subplot(3,2,6);
    plot(time,omega_z)
    hold on
    plot(time,filt_omega_z)
    legend('\omega_z','filtered \omega_z');
end

%% Calculation of angles using filtered data
% Angle calculated from accelerometer
roll_acc = atan2(-filt_a_y2,-filt_a_z2)*180/pi;
pitch_acc = atan2(filt_a_x2,sqrt(filt_a_y2.^2+filt_a_z2.^2))*180/pi;

% Angle initialization (first measure using accelerometer data only)
roll = zeros(1, length(time));
pitch = zeros(1, length(time));
roll(1) = roll_acc(1); % angle of rotation around x
pitch(1) = pitch_acc(1); % angle of rotation around y

% Complementary filtering using gyroscope and accelerometer data
coef1 = 0.98; % coef complemenary filter
coef2 = 1-coef1;
for t = 2:(length(time))
    roll(t) = coef1*(roll(t-1)+filt_omega_x(t)*(time(t)-time(t-1))) + coef2*roll_acc(t);
    pitch(t) = coef1*(pitch(t-1)+filt_omega_y(t)*(time(t)-time(t-1))) + coef2*pitch_acc(t);
    roll(t) = mod(roll(t)+180,360)-180;
    pitch(t) = mod(pitch(t)+180,360)-180;
end

if(plot_index == 1 || plot_index == 2)
    % plot accelerometer angles
    figure
    subplot(2,1,1)
    hold on
    plot(time, roll_acc,'r');
    plot(time, pitch_acc,'g');
    title('Tilt angles accelerometer'); xlabel('Time(s)'); ylabel('Angle(°)');
    legend('roll \phi','pitch \theta');
    hold off
    % plot data fusion angles (complementary filter)
    subplot(2,1,2)
    hold on
    plot(time, roll,'r');
    plot(time, pitch,'g');
    title('Tilt angles after complementary filter'); xlabel('Time(s)'); ylabel('Angle(°)');
    legend('roll \phi','pitch \theta');
    hold off
end

%% Remove gravity from accelerometer
% Calculate gravity in watch plan
g = [0; 0; -980];
new_g = zeros(3,length(time));
for t = 1:length(time)
    Rx = [1, 0, 0; 
          0, cosd(-roll_acc(t)), -sind(-roll_acc(t));
          0, sind(-roll_acc(t)), cosd(-roll_acc(t))];

    Ry = [cosd(-pitch_acc(t)), 0, sind(-pitch_acc(t));
          0, 1, 0;
          -sind(-pitch_acc(t)), 0, cosd(-pitch_acc(t))];
%     Rx = [1, 0, 0; 
%           0, cosd(roll(t)), -sind(roll(t));
%           0, sind(roll(t)), cosd(roll(t))];
% 
%     Ry = [cosd(pitch(t)), 0, sind(pitch(t));
%           0, 1, 0;
%           -sind(pitch(t)), 0, cosd(pitch(t))];
      
    new_g(:,t) = Rx*Ry*g;
end
a_x_gravity = filt_a_x - new_g(1,:).';
a_y_gravity = filt_a_y - new_g(2,:).';
a_z_gravity = filt_a_z - new_g(3,:).';

if(plot_index == 1 || plot_index == 2)
    % Plot data for general overview
    figure
    % Plot acceleromter data
    subplot(3,1,1);
    hold on
    plot(time, filt_a_x, 'r');
    plot(time, filt_a_y, 'g');
    plot(time, filt_a_z, 'b');
    v = axis;
    title('Accelerometer'); xlabel('time(s)'); ylabel('Acceleration(mg)');
    legend('x','y','z');
    hold off
    % Plot acceleromter without gravity data
    subplot(3,1,2);
    hold on
    plot(time, a_x_gravity, 'r');
    plot(time, a_y_gravity, 'g');
    plot(time, a_z_gravity, 'b');
    axis(v)
    title('Accelerometer with removed gravity'); xlabel('time(s)'); ylabel('Acceleration(mg)');
    legend('x','y','z');
    hold off
    % Plot gyroscope data
    subplot(3,1,3);
    hold on
    plot(time, omega_x, 'r');
    plot(time, omega_y, 'g');
    plot(time, omega_z, 'b');
    title('Gyroscope'); xlabel('Time(s)'); ylabel('Angular rate(°/s)');
    legend('x','y','z');
    hold off
end
%% Analyse datas 
% put here start and end time for data analysis (in seconds)
t_start = time(1);% t_start = time(1);
t_end = time(end);% t_end = time(end);
% index for start and end time
ts = 0; te = 0;
for t = 1:length(time)
    if(time(t)>t_start && ts==0)
        ts = t;
    end
    if(time(t)>=t_end && te==0)
        te = t-1;
    end
end
time(ts)
time(te)

% Acceleration maximums
max_xy_acc = [max(sqrt(a_x_gravity(ts:te).^2+a_y_gravity(ts:te).^2)),...
              mean(sqrt(a_x_gravity(ts:te).^2+a_y_gravity(ts:te).^2)),...
              max(abs(filt_omega_z(ts:te))),...
              mean(abs(filt_omega_z(ts:te)))];
          
max_min_data_gravity = [max_xy_acc,...
                max(a_x_gravity(ts:te)),    min(a_x_gravity(ts:te)),...
                max(a_y_gravity(ts:te)),    min(a_y_gravity(ts:te)),...
                max(a_z_gravity(ts:te)),    min(a_z_gravity(ts:te)),...
                max(filt_omega_x(ts:te)),min(filt_omega_x(ts:te)),...
                max(filt_omega_y(ts:te)),min(filt_omega_y(ts:te)),...
                max(filt_omega_z(ts:te)),min(filt_omega_z(ts:te))];

% Angles mean and variation 
figure
boxplot([roll_acc(ts:te), pitch_acc(ts:te)],{'roll','pitch'});
