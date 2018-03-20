%% DITHER
%  Juri Lukkarila 2017

%% 1kHz Sine 0 dB with dither

% read audiofile
[x, Fs] = audioread('1kHz-0-dither.wav');
info    = audioinfo('1kHz-0-dither.wav');

N = max(size(x));               % samples
Ts = 1/Fs;                      % sample time 
t = 0:Ts:(N-1)*Ts;              % time vector

%% waveform
figure();
plot(t, x); grid on;
title({'16-bit 1 kHz sine 0 dB with dither';'Waveform'});
xlabel('Time (ms)'); 
axis([0 0.01 -1.1 1.1]);
set(gca,'XTick', 0:0.001:0.01);
set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-0db', '-dpng', '-r300'); 

%% stem
figure();
stem(t(1:0.002*Fs), x(1:0.002*Fs), 'filled','MarkerSize', 5); grid on;
title({'16-bit 1 kHz sine 0 dB with dither';'Sample values'});
xlabel('Time (ms)');
axis([0 0.002 -1.1 1.1]);
set(gca,'XTick', 0:0.0002:0.002);
set(gca,'XTickLabel', 0:0.2:2);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz_0db_stem', '-dpng', '-r300'); 

%% level
l = 0.1;
figure();
plot(t(1:l*Fs), 20*log10(abs(x(1:l*Fs)))); grid on;
title({'16-bit 1 kHz sine 0 dB with dither';'Level (dBFS)'});
xlabel('Time (s)'); 
ylabel('Level (dBFS)');
axis([0 l -10 3]);
set(gca,'XTick', 0:l/10:l);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-0db_level', '-dpng', '-r300'); 

%% FFT
Xf = 2*abs((fft(x)/N));         % real part
Xf = 20*log10(Xf);              % dB
f0 = Fs/N;                      % frequency resolution (Hz)
f = 0:f0:(N-1)*f0;              % frequency vector

% only draw to 24 kHz -> fft bin index for 24 kHz when Fs = 48 kHz:
index = max(size(x))/2

% spectrum
figure();
semilogx(f(1:720000), Xf(1:720000)); grid on;
title({'16-bit 1 kHz sine 0 dB with dither';'Frequency spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([20 24000 -170 3]);
set(gca,'XTick',     [20 30 60 125 250 500 1000 2000 4000 8000 16000 24000])
set(gca,'XTickLabel',{20 30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-0db_spectrum', '-dpng', '-r300'); 

%% 1kHz Sine -105 dB with dither

% read audiofile
[x, Fs] = audioread('1kHz-105-dither.wav');
info2   = audioinfo('1kHz-105-dither.wav');

N = max(size(x));               % samples
Ts = 1/Fs;                      % sample time  
t = 0:Ts:(N-1)*Ts;              % time vector

%% waveform
figure();
plot(t, x); grid on;
title({'16-bit 1 kHz sine -105 dB with dither';'Waveform'});
xlabel('Time (ms)');
axis([0 0.01 -0.0005 0.0005]);
set(gca,'XTick', 0:0.001:0.01);
set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-105db', '-dpng', '-r300'); 

%% stem plot
figure();
stem(t(1:0.002*Fs), x(1:0.002*Fs), 'filled','MarkerSize', 5); grid on;
title({'16-bit 1 kHz sine -105 dB';'Sample values'});
xlabel('Time (ms)');
axis([0 0.002 -0.0005 0.0005]);
set(gca,'XTick', 0:0.0002:0.002);
set(gca,'XTickLabel', 0:0.2:2);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-105db_stem', '-dpng', '-r300'); 

%% level
figure();
plot(t, 20*log10(abs(x))); grid on;
title({'16-bit 1 kHz sine -105 dB with dither';'Level (dBFS)'});
xlabel('Time (s)'); 
ylabel('Level (dBFS)');
axis([0 1 -90 3]);
set(gca,'XTick', 0:0.1:1);
%set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-105db_level', '-dpng', '-r300'); 

%% normalize
xmin = min(x);
xmax = max(x);
if abs(xmin) >= xmax
    xnorm = x./abs(xmin);
else
    xnorm = x./abs(xmax);
end

% export
audiowrite('1kHz_normalized.wav', xnorm, Fs);

%% FFT
Xf = 2*abs((fft(x)/N));         % real part
Xf = 20*log10(Xf);              % dB
f0 = Fs/N;                      % frequency resolution (Hz)
f = 0:f0:(N-1)*f0;              % frequency vector

% spectrum
figure();
semilogx(f(1:720000), Xf(1:720000)); grid on;
title({'16-bit 1 kHz sine -105 dB with dither';'Frequency spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([20 24000 -170 3]);
set(gca,'XTick',     [20 30 60 125 250 500 1000 2000 4000 8000 16000 24000])
set(gca,'XTickLabel',{20 30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-105db_spectrum', '-dpng', '-r300'); 

%% spectrum zoom
figure();
semilogx(f(1:720000), Xf(1:720000)); grid on;
title({'16-bit 1 kHz sine -105 dB with dither';'Frequency spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([20 24000 -170 -100]);
set(gca,'XTick',     [20 30 60 125 250 500 1000 2000 4000 8000 16000 24000])
set(gca,'XTickLabel',{20 30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, '1kHz-105db_spectrum2', '-dpng', '-r300'); 

%% spectrogram
window = round(Fs/20);      % divide to approx. 100 windows
if mod(window, 2) ~= 0      % if odd
    window = window + 1;    
end
overlap = window/2;         % window overlap
nfft    = 2^15;             % fft points

figure();
spectrogram(x(1:2*Fs), hamming(window), overlap, nfft, Fs,...
            'yaxis','MinThreshold', -160, 'power');
title({'16-bit 1 kHz sine -105 dB with dither';'Spectrogram'});
set(gca,'YScale','log');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
axis([1 2 0.02 24]); 
set(gca,'YTick',[0.02 0.03 0.060 0.125 0.25 0.5 1 2 4 8 16 24]); 
set(gca,'YTickLabel',{20 30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
c = colorbar; c.Label.String = 'Magnitude (dB)';
view([75 50]);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 16 9])
print(gcf,'1kHz-105db_spectrogram_3D', '-dpng', '-r300'); 

%% Dither strong

% read audiofile
[x, Fs] = audioread('dither_16_strong.wav');
info4   = audioinfo('dither_16_strong.wav');

N = max(size(x));           % samples
Ts = 1/Fs;                  % sample time 
t = 0:Ts:(N-1)*Ts;          % time vector
 

%% waveform
figure();
plot(t, x); grid on;
title({'16-bit dither';'Waveform'});
xlabel('Time (ms)');
axis([0 0.01 -0.05 0.05]);
set(gca,'XTick', 0:0.001:0.01);
set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'dither', '-dpng', '-r300'); 

%% waveform zoom
figure();
plot(t, x); grid on;
title({'16-bit dither';'Waveform'});
xlabel('Time (ms)');
axis([0 0.005 -0.0005 0.0005]);
set(gca,'XTick', 0:0.001:0.01);
set(gca,'XTickLabel',{0 1 2 3 4 5 6 7 8 9 10});
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'dither_zoom', '-dpng', '-r300');

%% stem
figure();
stem(t(1:0.002*Fs), x(1:0.002*Fs), 'filled','MarkerSize', 5); grid on;
title({'16-bit dither';'Sample values'});
xlabel('Time (ms)');
axis([0 0.002 -0.0005 0.0005]);
set(gca,'XTick', 0:0.0002:0.002);
set(gca,'XTickLabel', 0:0.2:2);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'dither_stem', '-dpng', '-r300');

%% level
figure();
plot(t, 20*log10(abs(x))); grid on;
title({'16-bit dither';'Level (dB)'});
xlabel('Time (s)'); 
ylabel('Level (dB)');
axis([0 1 -90 3]);
set(gca,'XTick', 0:0.1:1);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'dither_level', '-dpng', '-r300'); 

%% normalize
xmin = min(x);
xmax = max(x);
if abs(xmin) >= xmax
    xnorm = x./abs(xmin);
else
    xnorm = x./abs(xmax);
end

% export
audiowrite('dither_normalized.wav', xnorm, Fs);

%% FFT
Xf = 2*abs((fft(x)/N));         % real part
Xf = 20*log10(Xf);              % dB
f0 = Fs/N;                      % frequency resolution (Hz)
f = 0:f0:(N-1)*f0;              % frequency vector

% spectrum
figure();
semilogx(f(1:720000), Xf(1:720000)); grid on;
title({'16-bit dither, strong setting';'Frequency spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([20 24000 -160 3]);
set(gca,'XTick',     [30 60 125 250 500 1000 2000 4000 8000 16000])
set(gca,'XTickLabel',{30 60 125 250 500 '1k' '2k' '4k' '8k' '16k'})
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
print(gcf, 'dither_spectrum', '-dpng', '-r300'); 