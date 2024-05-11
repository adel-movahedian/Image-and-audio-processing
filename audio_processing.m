%--------------------------------------------------------------------------
% 
%                               -Q1-
%
% --------------------------------------------------------------------------
%% Q 1 - discrete
clear;clc;
n = -10:10 ;
x_1 = sin (pi*n/2) + cos( pi * n) ;
subplot(2,2,1);
stem(n,x_1,'filled');
ylim([-3,3]);
grid minor;
title('$x_1[n]$','Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');
x_2 = zeros(1,length(n));
for k = -5:5
    x_2 = x_2 + dirac(n - 2*k);
end
idx = find(x_2 == Inf);
x_2(idx) = 1;
subplot(2,2,2);
stem(n,x_2,'filled','color' , 'r');
ylim([-1,2]);
grid minor;
title('$x_2[n]$','Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');
x_3 = exp(-(n.^2)/18);
subplot(2,2,3);
stem(n,x_3,'filled','color' , 'g');
ylim([-1,2]);
grid minor;
title('$x_3[n]$','Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');
N = 4;
x_4 = (mod(n,N) <(N/2));
subplot(2,2,4);
stem(n,x_4,'filled','color' , 'b');
ylim([-1,2]);
grid minor;
title('$x_4[n]$','Interpreter','latex');
xlabel('$n$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');

%% Q 1 - continuous
dt = 0.01;
T = 1;
t = -10:dt:10;
x_1 = sinc(t-pi);

subplot(3,2,1);
plot(t,x_1,'Color','b');
ylim([-0.5,1.5]);
grid minor;
title('$x_1(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');

subplot(3,2,2);
x2=0;
for j=1 : 19
x2 = x2 + (t-(j-10)).*(t-(j-10)>=0) - (t-0.5-(j-10)).*(t-0.5-(j-10)>=0) - 0.5*heaviside(t-0.5-j+10);
end
plot(t,x2,'Color','g')
ylim([-0.2,0.7]);
grid minor;
title('$x_2(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');

x_3 = dirac(t -10*T) ;
for n = -5:4
    x_3 = x_3 + dirac(t -2*n*T) - dirac(t -(2*n+1)*T);
end
index_negative = find(x_3 == -Inf);
x_3(index_negative) =-1;
index_positive = find(x_3 == Inf);
x_3(index_positive) = 1;
subplot(3,2,3)
plot(t,x_3,'Color','r');
ylim([-1.2,1.2]);
grid minor;
title('$x_3(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');

subplot(3,2,4);
x4 = @(t) (abs(mod(t,1))<1/4) + (abs(mod(t,1))>3/4) ;
plot(t,x4(t),'Color','m')
ylim([-0.2,1.2]);
grid minor;
title('$x_4(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');

x_5 = zeros(1,length(t));
for ta = 2:length(t)
   x_5(ta) = x_5(ta-1) + x_3(ta);
end
subplot(3,2,5)
plot(t,x_5,'-b');
ylim([-1.2,0.2]);
grid minor;
title('$x_5(t)$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$Amplitude$','Interpreter','latex');
%--------------------------------------------------------------------------
% 
%                               -Q2-
%
% --------------------------------------------------------------------------
%% Q 2-1 
clc; clear; close all;
[x, Fs] = audioread('Music.wav');
t = (0:length(x) - 1) / Fs;
plot(t, x);
title('Waveform');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Spectrogram');
soundsc(x,Fs);
%% Q 2-2
clc; clear; close all;
[x, Fs] = audioread('Music.wav');
n_d = round(0.2 * Fs);
a = 0.2;
x1 = (x(:, 1) + x(:, 2)) / 2;
y2 = zeros(length(x1) + n_d, 1);
y2(n_d + 1:end) = x1(1:end); 
y2 = a * y2;
y1 = x1;
y1(length(x1) + n_d) = 0;
y = y1 + y2;
soundsc(y, Fs);
t = (0:length(y) - 1) / Fs;
plot(t, y);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('echo Audio Signal');

%% Q 2-3 - uniform

clc; clear; close all;
[x, Fs] = audioread('Music.wav');
x1 = (x(:, 1) + x(:, 2)) / 2;
n1 = 0.2 * rand(length(x1), 1) - 0.1;
t = (0:length(x1) - 1) / Fs;
subplot(2, 1, 1);
plot(t, n1);
title('Noise');
ylim([-1,1]);
xlabel('Time (seconds)');
ylabel('Amplitude');
subplot(2, 1, 2);
plot(t, x1);
title('noised Audio Signal');
xlabel('Time (seconds)');
ylabel('Amplitude');
y = x1 + n1;
soundsc(y, Fs);

%% Q 2-3 - normal

clc; clear; close all;
[x, Fs] = audioread('Music.wav');
x1 = (x(:, 1) + x(:, 2)) / 2;
n2 = 0.01 * randn(length(x1), 1);
t = (0:length(x1) - 1) / Fs;
subplot(2, 1, 1);
plot(t, n2);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Generated Noise');
y = x1 + n2;
soundsc(y, Fs);
subplot(2, 1, 2);
plot(t, y);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('noised Audio Signal');

%% Q 2-4 
clc; clear; close all;[x, Fs] = audioread('Music.wav');
x1 = (x(:, 1) + x(:, 2)) / 2;
t = (1:length(x1) - 1) / Fs;
f = 2000 : (1000 / (length(x1) - 1)) : 3000;
Sin = 0.1 * sin(2 * pi * f);
y = x1' + Sin;
soundsc(y, Fs);
audiowrite('MusicPlusSin.wav', y, Fs);

%--------------------------------------------------------------------------
% 
%                               -Q3-
%
% --------------------------------------------------------------------------
%% Q 3 -1
n=0:5;
n2=-2:2;
x=[2 4 6 -2 1 2];
l=[1/3 2/3 1 2/3 1/3];
figure
subplot(2,2,1);
stem(n,x);
title('$x[n]$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-3 7]);
xlim([-1 6]);
grid minor;
hold on;
plot(n,x,"--",LineWidth=1,Color="r")
x_3 = zeros(1, 3*length(x)- 2);
x_3(1:3:end)= x;
x_3(2:3:end) = 0;
x_3(3:3:end) = 0;
subplot(2,2,2);
stem(x_3);
title('$x_3[n]$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-3 7]);
xlim([-1 18]);
grid minor;
hold on;
linex_3 = zeros(1, 3*length(x)- 2);
linex_3(1:3:end)= x;
linex_3(2:3:end) = (2/3)*x(1:end-1)+(1/3)*x(2:end);
linex_3(3:3:end) = (1/3)*x(1:end-1)+(2/3)*x(2:end);
x_points = [1:length(x_3)/16:16];
y_points = linex_3;
plot(x_points,y_points,"--",LineWidth=1,Color="r")
subplot(2,2,3);
stem(n2,l);
title('$l[n]$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-0.5,2]);
xlim([-3 3]);
grid minor;
hold on;
x_points = [-2:length(l)/5:2];
y_points = l;
plot(x_points,y_points,"--",LineWidth=1,Color="r")
lin_intp = conv(x_3, l);
subplot(2,2,4);
n2 =1:20;
stem(n2, lin_intp);
title('$convolved\ signal$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-3,7]);
xlim([0 20]);
grid minor;
hold on;
x_points = n2;
y_points = lin_intp;
plot(x_points,y_points,"--",LineWidth=1,Color="r")
%% Q 3 -2
n=0:5;
n2=-2:2;
x=[2 4 6 -2 1 2];
l=[1/3 2/3 1 2/3 1/3];
figure
subplot(2,2,1);
stem(n,x);
title('$x[n]$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-3 7]);
xlim([-1 6]);
grid minor;
hold on;
plot(n,x,"--",LineWidth=1,Color="r")
x_2 = zeros(1, 2*length(x)- 1);
x_2(1:2:end)= x;
x_2(2:2:end) = 0;
subplot(2,2,2);
stem(x_2);
title('$x_2[n]$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-3 7]);
xlim([-1 14]);
grid minor;
hold on;
linex_2 = zeros(1, 2*length(x)- 1);
linex_2(1:2:end)= x;
linex_2(2:2:end) = (1/2)*x(1:end-1)+(1/2)*x(2:end);
x_points = [1:length(x_2)/11:11];
y_points = linex_2;
plot(x_points,y_points,"--",LineWidth=1,Color="r")
subplot(2,2,3);
stem(n2,l);
title('$l[n]$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-0.5,2]);
xlim([-3 3]);
grid minor;
hold on;
x_points = [-2:length(l)/5:2];
y_points = l;
plot(x_points,y_points,"--",LineWidth=1,Color="r")
lin_intp = conv(x_2, l);
subplot(2,2,4);
n2 =1:15;
stem(n2, lin_intp);
title('$convolved\ signal$','Interpreter','latex');
xlabel('n');
ylabel('Amplitude');
ylim([-3,8]);
xlim([0 16]);
grid minor;
hold on;
x_points = n2;
y_points = lin_intp;
plot(x_points,y_points,"--",LineWidth=1,Color="r")
