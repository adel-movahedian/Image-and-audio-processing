%% Question 1
clear; clc; close all;
fs = 200;
t = 0:1/fs:5; 
A_values = arrayfun(@A, t);
n_values = n(t);
f_values = A_values + n_values;
figure;
subplot(2, 1, 1);
plot(t, A_values, 'r', t, f_values, 'g');
xlabel('Time');
ylabel('Amplitude');
legend('A(t)', 'f(t)');
title('A(t) & f(t)');
grid on;
%-------
f_fft = fft(f_values);
f_shifted = fftshift(f_fft);
frequencies = linspace(-fs/2, fs/2, length(t));
subplot(2, 1, 2);
plot(frequencies, abs(f_shifted), 'g');
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('F(f)');
grid on;
%-------
fc = 50;
H_idl = (abs(frequencies) <= fc);
figure;
plot(frequencies, H_idl);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('low pass ideal filter H(f)');
xlim([-100,100]);
ylim([-0.5,1.5]);
grid on;
%-------
Y_fft = f_shifted .* H_idl;
figure;
plot(frequencies, abs(Y_fft), 'r');
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Y(f)');
grid on;
%-------
n_shifted = fftshift(fft(n_values));
n_shifted_withoutcos = n_shifted .* (1 - H_idl);  
Y_shiftedInverse = ifftshift(Y_fft);
n_shifted_withoutcosInverse = ifftshift(n_shifted_withoutcos);
A_time = ifft(Y_shiftedInverse);
n_time_withoutcos = ifft(n_shifted_withoutcosInverse);
output_time = n_time_withoutcos + A_time - n(t) ; 
figure;
plot(t, real(output_time), 'r'); 
xlabel('Time');
ylabel('Amplitude');
title('regenerated A');
grid on;
%% bonus question 1
A_values = arrayfun(@A, t);
n_values = zeros(size(t));
n_values(t > 0 & t < 5) = 3 * exp(-10 * t(t > 0 & t < 5));
f_values = A_values + n_values;
figure;
subplot(2, 1, 1);
plot(t, A_values, 'b', t, f_values, 'r');
xlabel('Time');
ylabel('Amplitude');
legend('A(t)', 'f(t)');
title('A(t) & f(t)');
grid on;
%-------
f_fft = fft(f_values);
f_shifted = fftshift(f_fft);
frequencies = linspace(-fs/2, fs/2, length(t));
subplot(2, 1, 2);
plot(frequencies, abs(f_shifted), 'b');
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('F(f)');
grid on;
%-------
fc = 50;
H_idl = (abs(frequencies) <= fc);
figure;
plot(frequencies, H_idl);
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('low pass ideal filter H(f)');
ylim([-1,2]);
grid on;
%-------
Y_fft = f_shifted .* H_idl;
figure;
plot(frequencies, abs(Y_fft), 'r');
xlabel('Frequency (Hz)');
ylabel('magnitude');
title('Y(f)');
grid on;
%-------
Y_fft = f_shifted .* H_idl;
A_time = ifft(ifftshift(Y_fft)) - n_values ;
figure;
plot(t, real(A_time), 'g');
xlabel('Time');
ylabel('Amplitude');
title('regenerated A');
grid on;
%% Question 2
clear; clc; close all;
fs = 300;  
Ts = 1/fs;
t = 0:Ts:5; 
A1_t = signal1(t);
A2_t = signal2(t);
w1_t = A1_t .* cos(2*pi*30*t);
w2_t = A2_t .* cos(2*pi*60*t);
w_t =  w1_t + w2_t ;
mod_signal1 = 2 * cos(2*pi*30*t) .* w_t;
mod_signal2 = 2 * cos(2*pi*60*t) .* w_t;
mod_signal_f_1 = fft(mod_signal1);
mod_signal_f_1 = fftshift(mod_signal_f_1);
mod_signal_f_2 = fft(mod_signal2);
mod_signal_f_2 = fftshift(mod_signal_f_2);
%------
frequencies = linspace(-fs/2, fs/2, length(t));
fc = 10;
H_idl = (abs(frequencies) <= fc);
filtered_signal_f_1 = mod_signal_f_1 .* H_idl;
demod_signal_t_1 = ifft(ifftshift(filtered_signal_f_1));
filtered_signal_f_2 = mod_signal_f_2 .* H_idl;
demod_signal_t_2 = ifft(ifftshift(filtered_signal_f_2));
%------
figure;
plot(t, A1_t + A2_t, 'b');
title(' A1(t)+A2(t)');
xlabel('Time');
ylim([-1,3]);
grid on
%--------
figure;
plot(t, A1_t, 'b', t, A2_t, 'r');
title('Original  A1(t) & A2(t)');
xlabel('Time');
ylim([-1,2]);
grid on
%-------
figure;
plot(t, w_t);
title(' w(t)');
xlabel('Time');
ylim([-1.5,2.5]);
grid on
%-------
figure;
plot(t , mod_signal1 , 'r' );
title('modulated Signal 1');
xlabel('time');
ylim([-2.5,4.5]);
grid on
%------
figure;
plot(frequencies, abs(mod_signal_f_1) , 'black');
title('modulated Signal 1 (frequency domain)');
xlabel('frequency (Hz)');
grid on
%------
figure;
plot(frequencies , abs(filtered_signal_f_1),'g');
title('filtered Signal 1 (frequency domain)');
xlabel('frequency (Hz)');
grid on
%------
figure;
plot(t , real(demod_signal_t_1),'r');
title('demodulated signal A1(T)');
xlabel('time');
grid on
%------
figure;
plot(t , mod_signal2 );
title('modulated Signal 2');
xlabel('time');
ylim([-2.5,4.5]);
grid on
%------
figure;
plot(frequencies, abs(mod_signal_f_2) , 'm');
title('modulated Signal 2 (frequency domain)');
xlabel('frequency (Hz)');
grid on
%------
figure;
plot(frequencies , abs(filtered_signal_f_2),'r');
title('filtered Signal 2 (frequency domain)');
xlabel('frequency (Hz)');
grid on
%------
figure;
plot(t , real(demod_signal_t_2),'b');
title('demodulated signal A2(T)');
xlabel('time');
grid on

%% Functions
% signal1(t)
function y = signal1(t)
    y = zeros(size(t));
    y((t > 1 & t < 2) | (t > 3 & t < 4)) = 1;
end
% signal2(t)
function y = signal2(t)
    y = zeros(size(t));
    y((t > 0.5 & t < 1.5) | (t > 2.5 & t < 3.5)) = 1;
end
% A(t)
function y = A(t)
    if t <= 2
        y = 1;
    elseif t > 2
        y = 0.1;
    else
        y = 0;
    end
end
% n(t)
function y = n(t)
    y = zeros(size(t));
    y(t > 0 & t < 5) = 0.5 * cos(2 * pi * 50 * t(t > 0 & t < 5));
end
