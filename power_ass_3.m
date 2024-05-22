close all;
clear all;
clc;
pair_poles = 8/2;
S_min = -1;
S_max = 2;
R1 = 0.085;
R2_prime = 0.067;
X1 = 0.196i;
X2_prime = 0.161i;
X_muo = 6.65i;
Vph = 220-127.0170592i;
% a) at Rated conditions
F = 50; %Hz
Vth = Vph*(X_muo/(X_muo+X1+R1));
Zth = (X_muo*(X1+R1))/(X1+X_muo+R1); %ohm
ns = (60*F)/pair_poles; %rpm
Ws = (2*pi*ns)/60; %rad/s
n = [(ns)*(1-S_max):0.1:(ns)*(1-S_min)]; %rpm
%caculating the torque
Torque = zeros(1,length(n));
for i=1:length(n)
    Torque(i) = (3*(abs((Vth))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth+(0.067*(ns/(ns-n(i))))+0.161i))^2);
end
%plotting the torque Vs speed
figure('Name','The torque-speed characteristics for a slip range from 2 to-1 at rated conditions.');
plot(n,Torque,'r','linewidth',1);
xlim([-1000 1600]);
title("Torque Vs Speed plot for rated conditions");
xlabel("Speed (rpm)");
ylabel("Torque (N.m)");
grid on;

%b)at percents of Vrated(same parameters as required (a)) but:
%only multiply Vth by percents required

%caculating torques
Torque_80percent_Vrated = zeros(1,length(n));
Torque_60percent_Vrated = zeros(1,length(n));
Torque_40percent_Vrated = zeros(1,length(n));
Torque_20percent_Vrated = zeros(1,length(n));

for i=1:length(n)
    Vth8 = 0.8*Vth;
    Torque_80percent_Vrated(i) = (3*(abs((Vth8))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth+(0.067*(ns/(ns-n(i))))+0.161i))^2);
    Vth6 = 0.6*Vth;
    Torque_60percent_Vrated(i) = (3*(abs((Vth6))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth+(0.067*(ns/(ns-n(i))))+0.161i))^2);
    Vth4 = 0.4*Vth;
    Torque_40percent_Vrated(i) = (3*(abs((Vth4))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth+(0.067*(ns/(ns-n(i))))+0.161i))^2);
    Vth2 = 0.2*Vth;
    Torque_20percent_Vrated(i) = (3*(abs((Vth2))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth+(0.067*(ns/(ns-n(i))))+0.161i))^2);
end

% plotting torques Vs speed
figure('Name' , 'The torque-speed characteristics for a slip range from 2 to -1 at percents of Vrated');
plot(n,Torque_80percent_Vrated,'r','linewidth',1);
hold on;
plot(n,Torque_60percent_Vrated,'g','linewidth',1);
hold on;
plot(n,Torque_40percent_Vrated,'y','linewidth',1);
hold on;
plot(n,Torque_20percent_Vrated,'k','linewidth',1);
hold off;
title('Torque Vs Speed plot for percent of rated voltage');
xlabel("Speed (rpm)");
ylabel("Torque (N.m)");
grid on;
legend('V = 0.8*Vrated','V = 0.6*Vrated','V = 0.4*Vrated','V = 0.2*Vrated','Location','NorthEast');
xlim([-1000 1600]);
grid on;

%c)at different frequencies
%calculating the changing parameters
F15 = 15;
F25 = 25;
F65 = 65;

X1_15Hz = X1*(F15/F);
X2_prime_15Hz = X2_prime*(F15/F);
X_muo_15Hz = X_muo*(F15/F);

X1_25Hz = X1*(F25/F);
X2_prime_25Hz = X2_prime*(F25/F);
X_muo_25Hz = X_muo*(F25/F);

X1_65Hz = X1*(F65/F);
X2_prime_65Hz = X2_prime*(F65/F);
X_muo_65Hz = X_muo*(F65/F);

Zth = (X_muo*(X1+R1))/(X1+X_muo+R1);

Vth_15Hz = Vph*(X_muo_15Hz/(X_muo_15Hz+X1_15Hz+R1));
Vth_25Hz = Vph*(X_muo_25Hz/(X_muo_25Hz+X1_25Hz+R1));
Vth_65Hz = Vph*(X_muo_65Hz/(X_muo_65Hz+X1_65Hz+R1));

Zth_15Hz = (X_muo_15Hz*(X1_15Hz+R1))/(X1_15Hz+X_muo_15Hz+R1);
Zth_25Hz = (X_muo_25Hz*(X1_25Hz+R1))/(X1_25Hz+X_muo_25Hz+R1);
Zth_65Hz = (X_muo_25Hz*(X1_25Hz+R1))/(X1_25Hz+X_muo_25Hz+R1);

ns_15Hz = (60*F15)/pair_poles;
ns_25Hz = (60*F25)/pair_poles;
ns_65Hz = (60*F65)/pair_poles;

Ws_15Hz = (2*pi*ns_15Hz)/60;
Ws_25Hz = (2*pi*ns_25Hz)/60;
Ws_65Hz = (2*pi*ns_65Hz)/60;

n_15HZ = [(ns_15Hz)*(1-S_max):0.1:(ns_15Hz)*(1-S_min)];
n_25HZ = [(ns_25Hz)*(1-S_max):0.1:(ns_25Hz)*(1-S_min)];
n_65HZ = [(ns_65Hz)*(1-S_max):0.1:(ns_65Hz)*(1-S_min)];
% calculating torques
Torque_15Hz = zeros(1,length(n_15HZ));
Torque_25Hz = zeros(1,length(n_25HZ));
Torque_65Hz = zeros(1,length(n_65HZ));

for i=1:length(n_15HZ)
    Torque_15Hz(i) = (3*(abs((Vth_15Hz))^2)*0.067*(ns_15Hz/(ns_15Hz-n_15HZ(i)))) / (Ws_15Hz*(abs(Zth_15Hz+(0.067*(ns_15Hz/(ns_15Hz-n_15HZ(i))))+0.161i))^2);
end

for i=1:length(n_25HZ)
    Torque_25Hz(i) = (3*(abs((Vth_25Hz))^2)*0.067*(ns_25Hz/(ns_25Hz-n_25HZ(i)))) / (Ws_25Hz*(abs(Zth_25Hz+(0.067*(ns_25Hz/(ns_25Hz-n_25HZ(i))))+0.161i))^2);
end

for i=1:length(n_65HZ)
    Torque_65Hz(i) = (3*(abs((Vth_65Hz))^2)*0.067*(ns_65Hz/(ns_65Hz-n_65HZ(i)))) / (Ws_65Hz*(abs(Zth_65Hz+(0.067*(ns_65Hz/(ns_65Hz-n_65HZ(i))))+0.161i))^2);
end

%plotting torques Vs speed
figure('Name','The torque-speed characteristics for a slip range from 2 to-1 at different frequencies.');
plot(n_15HZ,Torque_15Hz,'r','linewidth',1);
hold on;
plot(n_25HZ,Torque_25Hz,'g','linewidth',1);
hold on;
plot(n,Torque,'y','linewidth',1);
hold on;
plot(n_65HZ,Torque_65Hz,'k','linewidth',1);
hold off;
title("Torque Vs Speed plot for different frequencies");
xlabel("Speed (rpm)");
ylabel("Torque (N.m)");
legend('F = 15Hz','F = 25Hz','F = 50Hz','F = 65Hz','Location','NorthEast');
xlim([-1000 2000]);
grid on;

% d) at different values for R1
%calculating the changing parameters
Vth_2R1 = 216.631639-118.0011781i;
Vth_3R1 = 217.9951021-115.2599516i;
Vth_4R1 = 219.2881314-112.489845i;
Zth_2R1 = 0.1603063428+0.1943692782i;
Zth_3R1 = 0.240274428+0.1993382967i;
Zth_4R1 = 0.3200210475+0.2062820853i;
%calculating torques
Torque_2R1 = zeros(1,length(n));
Torque_3R1 = zeros(1,length(n));
Torque_4R1 = zeros(1,length(n));
for i=1:length(n)
    Torque_2R1(i) = (3*(abs((Vth_2R1))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth_2R1+(0.067*(ns/(ns-n(i))))+0.161i))^2);
    Torque_3R1(i) = (3*(abs((Vth_3R1))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth_3R1+(0.067*(ns/(ns-n(i))))+0.161i))^2);
    Torque_4R1(i) = (3*(abs((Vth_4R1))^2)*0.067*(ns/(ns-n(i)))) / (Ws*(abs(Zth_4R1+(0.067*(ns/(ns-n(i))))+0.161i))^2);
end
%plotting troques Vs speed
figure('Name','The torque-speed characteristics for a slip range from 2 to-1 at circuit parameteres.');
plot(n,Torque_2R1,'b','linewidth',1);
hold on;
plot(n,Torque_3R1,'g','linewidth',1);
hold on;
plot(n,Torque_4R1,'k','linewidth',1);
hold off;
title("Torque Vs Speed plot for different circuit parameters");
xlabel("Speed (rpm)");
ylabel("Torque (N.m)");
legend('R1 = 2R1','R1 = 3R1','R1 = 4R1','Location','NorthEast');
xlim([-1000 1600]);
grid on;