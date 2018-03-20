%% SERATO DJ CONTROL SIGNAL ANALYSIS
%  Juri Lukkarila 2016

clear all;

% read audio file, 10 sec sample without the silence in the beginning
[xt, Fs] = audioread('Serato_Control_CD.wav', [1363, 1363+10*44100]); 

N = max(size(xt));          % samples
Ts = 1/Fs;                  % sample time 
t = 0:Ts:(N-1)*Ts;          % time vector

xt = xt/max(xt(:,1));       % normalize to 1

xleft = xt(1:N)';           % left channel
xright = xt(1:N, 2);        % right channel

%% Time plots

figure();
plot(t, xleft); grid on;
title({'Serato Control Signal';'left channel'});
xlabel('Time (ms)');
axis([0 0.1 -1.05 1.05]);
set(gca,'XTick',0:0.01:0.1);
set(gca,'XTickLabel',{0 10 20 30 40 50 60 70 80 90 100});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_left', '-dpng', '-r300'); 

figure();
plot(t, xright, 'r');  grid on;
title({'Serato Control Signal';'right channel'});
xlabel('Time (ms)');
axis([0 0.1 -1.05 1.05]);
set(gca,'XTick',0:0.01:0.1);
set(gca,'XTickLabel',{0 10 20 30 40 50 60 70 80 90 100});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_right', '-dpng', '-r300'); 

figure();
plot(t, xt); grid on;
title({'Serato Control Signal';'stereo'});
xlabel('Time (ms)');
axis([0 0.1 -1.05 1.05]);
set(gca,'XTick',0:0.01:0.1);
set(gca,'XTickLabel',{0 10 20 30 40 50 60 70 80 90 100});
legend('Left', 'Right', 'Location', 'Best');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_stereo', '-dpng', '-r300'); 

figure();
plot(t, xt); grid on;
title({'Serato Control Signal';'stereo'});
xlabel('Time (ms)');
axis([0 0.01 -1.05 1.05]);
set(gca,'XTick', 0:0.001:0.01);
set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10});
legend('Left', 'Right', 'Location', 'Best');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_stereo2', '-dpng', '-r300'); 

%% Frequency plots

% FFT
Xf = 2*abs((fft(xleft)/N));     % real part
Xf = 20*log(Xf);                % dB

f0 = Fs/N;                      % frequency resolution (Hz)
f1 = 0:f0:(N-1)*f0;             % frequency vector

figure();
semilogx(f1, Xf); grid on;
title({'Serato Scratch Live Control Signal';'frequency spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([20 20000 -90 3]);
set(gca,'XTick',[20 50 100 200 500 1000 2000 5000 10000 20000])
set(gca,'XTickLabel',{20 50 100 200 500 '1k' '2k' '5k' '10k' '20k'})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_fft', '-dpng', '-r300'); 

figure(6);
plot(f1, Xf); grid on;
title({'Serato Scratch Live Control Signal';'frequency spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([900 1100 -90 3]);
set(gca,'XTick',[950 1000 1050])
set(gca,'XTickLabel',{950 1000 1050})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_fft2', '-dpng', '-r300'); 

%% Pure sine wave reference

% export 10 sec clip
audiowrite('Serato Control Signal clip.wav', xt, Fs);

% 1 kHz sine wave for reference
sinref = zeros(max(size(t)), 2);
sinref(:,1) = sin(2*pi*1000.*t);
sinref(:,2) = sin(2*pi*1000.*t);
audiowrite('1 kHz sine reference.wav', sinref, Fs)

%% Spectrograms

% reference sine
figure();
spectrogram(sinref(:,1),1024, 512,[],Fs,'yaxis','MinThreshold',-100)
title('1 kHz sine, 1024 sample window'); axis([0 10 0 2]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Sine_spectrogram1', '-dpng', '-r300'); 

figure();
spectrogram(sinref(:,1),2^12, 2^11,[],Fs,'yaxis','MinThreshold',-100)
title('1 kHz sine, 4096 sample window'); axis([0 10 0 2]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Sine_spectrogram2', '-dpng', '-r300'); 

% serato
figure();
spectrogram(xleft(1:441001),1024, 512,[],Fs,'yaxis', 'MinThreshold',-100)
title('Serato timecode, 1024 sample window'); axis([0 10 0 2]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_spectrogram1', '-dpng', '-r300'); 

figure();
spectrogram(xleft(1:441001),2^12, 2^11,[],Fs,'yaxis', 'MinThreshold',-100)
title('Serato timecode, 4096 sample window'); axis([0 2 0 2]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_spectrogram2', '-dpng', '-r300'); 

figure();
spectrogram(xleft,2^11, 2^10,[],Fs,'yaxis', 'MinThreshold',-100)
title('Serato timecode, 2048 sample window'); axis([0.01 0.2 0 2]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_spectrogram3', '-dpng', '-r300'); 

%% Phase shift

samples90deg = 44100/(4*1000); % samples for 1/4 wavelength for 1 kHz sine

% -> samples90deg = 11

% pad right channel 11 samples = 90 deg phase shift backwards
% -> waveforms should overlap
rightpad = padarray(xright,11,'pre'); % insert 11 zeros to beginning
padded = rightpad(1:441001);

figure(5);
plot(t, xleft, t, padded); grid on;
title({'Serato Scratch Live Control Signal';'padded right channel'});
xlabel('Time (ms)');
axis([0 0.02 -1.05 1.05]);
set(gca,'XTick', 0:0.001:0.02);
set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20});
legend('Left', 'Right', 'Location', 'Best');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_padded', '-dpng', '-r300'); 

%% Amplitude shft keying demodulation

figure();
plot(t, xleft); grid on;
title({'Serato Control Signal';'Amplitude shift keying'});
xlabel('Time (s)');
axis([0 0.2 -0.05 1.05]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_ASK', '-dpng', '-r300'); 

% 1 kHz period in samples
period = 44100/1000; % = 44.1

figure();
plot(xleft); grid on; 
title({'Serato Control Signal';'amplitude shift keying'});
xlabel('Samples');
axis([439 1101 -0.05 1.05]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_ASK2', '-dpng', '-r300'); 

% simple demodulation assuming every period encodes one bit
bit_array = []; a = 1;
for i = 2:max(size(xleft))
    % take one period of signal and convert to bit value
    if xleft(i-1) < 0 && xleft(i) >= 0
        period = xleft(a:i); a = i;         % one period
        if any(period > 0.8)                % has a value larger than 0.8
            bit_array = [bit_array 1];      % append 1
        else
            bit_array = [bit_array 0];      % append 0
        end
    end
end

% Convert bit array to square wave representation
square = ones(441,1)*bit_array;     % replicate each value 441 times
square = square(:);                 % array to column

figure();
s = 1:132000;                       % sample vector
plot(s, xleft(1:132000)); hold on
plot(s/10, square(1:132000), 'r', 'LineWidth', 0.7); grid on; 
title({'Serato Control Signal';'ASK demodulation'});
xlabel('Bits'); axis([1 2205 -0.05 1.05]);
set(gca,'XTick', 0:44.1:2205); set(gca,'XTickLabel',[]);
set(gca,'YTick', [0 1]); set(gca,'YTickLabel',[0 1]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_bitmask', '-dpng', '-r300'); 

%% Bit array

x = 0:1:64;
figure();
stem(x, bit_array(1:max(size(x))), 'filled'); grid on; 
title({'Serato Control Signal';'bit values'});
xlabel('bit'); axis([0 64 -0.1 1.1]);
set(gca,'XTick', 0:8:64);
set(gca,'YTick', [0 1]); set(gca,'YTickLabel',[0 1]);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_bitvalues', '-dpng', '-r300'); 

%% Direction

t = 0:1/60000:0.1;

right = sin(2*pi*1000.*t);
left = sin(2*pi*1000.*t - pi/2);

figure();
plot(t(16:121), left(16:121), 'b'); hold on;
plot(t, right, 'r'); grid on;
title({'Serato Control Signal';'direction forward'});
xlabel('Time (ms)');
axis([0 0.002 -1.05 1.05]);
set(gca,'XTick',0:0.00025:0.002);
set(gca,'XTickLabel',{0 0.25 0.5 0.75 1 1.25 1.5 1.75 2});
legend('Left', 'Right', 'Location', 'Best');
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'Serato_direction', '-dpng', '-r300'); 
