%% DITHER
%  Juri Lukkarila
%  2017

close all; clc; screen = get(0,'screensize');

%% 1kHz Sine 0 dB with dither

[x, Fs] = audioread('sine_1kHz_0db_dither.wav'); % read audiofile
info    = audioinfo('sine_1kHz_0db_dither.wav')  % print info

N  = max(size(x));                        % samples
Ts = 1/Fs;                                % sample time 
t  = 0:Ts:(N-1)*Ts;                       % time vector
ymin = min(x); yminlim = ymin + 0.05*ymin;
ymax = max(x); ymaxlim = ymax + 0.05*ymax;

% waveform
figure('Position', [0, 0, 0.75*screen(3), 0.75*screen(4)]);
subplot(2,3,1); 
l  = 0.01;  % length to plot
p  = 4;     % oversample factor
xr = resample(x,p,1);
tr = 0:(1/(p*Fs)):(p*N-1)*(1/(p*Fs));
plot(tr, xr); grid on;
axis([0 l yminlim ymaxlim]);
xlabel('Time (ms)');
ylabel('Amplitude');
set(gca,'XTick', 0:l/10:l);
set(gca,'XTickLabel',0:(l/10*1000):l*1000);
title({'16-bit 1 kHz sine 0 dB peak with dither';'Waveform'});
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db', '-dpng', '-r300'); 

% stem
subplot(2,3,2);
l = 0.002;  % length to plot
stem(t(1:l*Fs), x(1:l*Fs), 'filled','MarkerSize', 3); grid on;
axis([0 l yminlim ymaxlim]);
xlabel('Time (ms)');
ylabel('Amplitude');
set(gca,'XTick', 0:l/10:l);
set(gca,'XTickLabel',0:(l/10*1000):l*1000);
title({'16-bit 1 kHz sine 0 dB peak with dither';'Sample values'});
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db_stem', '-dpng', '-r300'); 

% level
l = 0.02;
rms = sqrt(sum(x.^2)/length(x)); rms_db = 20*log10(rms);
subplot(2,3,3);
plot(t(1:l*Fs), 20*log10(abs(x(1:l*Fs)))); grid on; hold on;
line([0 l],[rms_db rms_db],'Color','Red','LineWidth',1);
axis([0 l -10 3]);
xlabel('Time (ms)'); 
ylabel('Level (dBFS)');
set(gca,'XTick', 0:l/10:l);
set(gca,'XTickLabel',0:(l/10*1000):l*1000);
title({'16-bit 1 kHz sine 0 dB peak with dither';'Level'});
legend('Peak','RMS','Location','Best');
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db_level', '-dpng', '-r300'); 

% spectrum
Xf = 2*abs((fft(x)/N));         % FFT real part
Xf = 20*log10(Xf);              % dB
f0 = Fs/N;                      % frequency resolution (Hz)
f = 0:f0:(N-1)*f0;              % frequency vector
index = find(f >= Fs/2, 1);     % only draw to Fs/2
subplot(2,3,[4,6]);
semilogx(f(1:index), Xf(1:index)); grid on;
title({'16-bit 1 kHz sine 0 dB peak with dither';'Magnitude spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([30 24000 -180 10]);
set(gca,'XTick',     [30 60 125 250 500 1000 2000 4000 8000 16000 24000])
set(gca,'XTickLabel',{30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db_spectrum', '-dpng', '-r300');
print(gcf, './figures/dither_1kHz_0db_combined', '-dpng', '-r150');

%% 1kHz Sine -105 dB with dither

% read audiofile
[x, Fs] = audioread('sine_1kHz_105db_dither.wav');
info    = audioinfo('sine_1kHz_105db_dither.wav')

N = max(size(x));               % samples
Ts = 1/Fs;                      % sample time  
t = 0:Ts:(N-1)*Ts;              % time vector
ymin = min(x); yminlim = ymin + 0.05*ymin;
ymax = max(x); ymaxlim = ymax + 0.05*ymax;

% waveform
figure('Position', [0, 0, 0.75*screen(3), 0.75*screen(4)]);
subplot(2,3,1); 
l  = 0.01;  % length to plot
p  = 4;     % oversample factor
xr = resample(x,p,1);
tr = 0:(1/(p*Fs)):(p*N-1)*(1/(p*Fs));
plot(tr, xr); grid on;
axis([0 l yminlim ymaxlim]);
xlabel('Time (ms)');
ylabel('Amplitude');
set(gca,'XTick', 0:l/10:l);
set(gca,'XTickLabel',0:(l/10*1000):l*1000);
title({'16-bit 1 kHz sine -105 dB peak with dither';'Waveform'});
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db', '-dpng', '-r300'); 

% stem
subplot(2,3,2);
l = 0.002;  % length to plot
stem(t(1:l*Fs), x(1:l*Fs), 'filled','MarkerSize', 3); grid on;
axis([0 l yminlim ymaxlim]);
xlabel('Time (ms)');
ylabel('Amplitude');
set(gca,'XTick', 0:l/10:l);
set(gca,'XTickLabel',0:(l/10*1000):l*1000);
title({'16-bit 1 kHz sine -105 dB peak with dither';'Sample values'});
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db_stem', '-dpng', '-r300'); 

% level
l = 0.5;
rms = sqrt(sum(x.^2)/length(x)); rms_db = 20*log10(rms);
subplot(2,3,3);
plot(t(1:l*Fs), 20*log10(abs(x(1:l*Fs)))); grid on; hold on;
line([0 l],[rms_db rms_db],'Color','Red','LineWidth',1);
axis([0 l -80 -60]);
xlabel('Time (ms)'); 
ylabel('Level (dBFS)');
set(gca,'XTick', 0:l/10:l);
set(gca,'XTickLabel',0:(l/10*1000):l*1000);
title({'16-bit 1 kHz sine -105 dB peak with dither';'Level'});
legend('Peak','RMS','Location','Best');
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db_level', '-dpng', '-r300'); 

% spectrum
Xf = 2*abs((fft(x)/N));         % FFT real part
Xf = 20*log10(Xf);              % dB
f0 = Fs/N;                      % frequency resolution (Hz)
f = 0:f0:(N-1)*f0;              % frequency vector
index = find(f >= Fs/2, 1);     % only draw to Fs/2
subplot(2,3,[4,6]);
semilogx(f(1:index), Xf(1:index)); grid on;
title({'16-bit 1 kHz sine -105 dB peak with dither';'Magnitude spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([30 24000 -180 10]);
set(gca,'XTick',     [30 60 125 250 500 1000 2000 4000 8000 16000 24000])
set(gca,'XTickLabel',{30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
%set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 16 9])
%print(gcf, './figures/dither_16bit_1kHz_0db_spectrum', '-dpng', '-r300');
print(gcf,'./figures/dither_1kHz_105db_combined','-dpng','-r150');

%% Spectrum zoom

figure('Position', [0, 0, 0.75*screen(3), 0.75*screen(4)]);
semilogx(f(1:720000), Xf(1:720000)); grid on;
title({'16-bit 1 kHz sine -105 dB peak with dither';'Magnitude spectrum'});
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
axis([30 24000 -170 -80]);
set(gca,'XTick',     [30 60 125 250 500 1000 2000 4000 8000 16000 24000])
set(gca,'XTickLabel',{30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
print(gcf,'./figures/dither_16bit_1kHz_105db_spectrum_zoom','-dpng','-r150'); 

%% Spectrogram

window = round(Fs/20);        % divide to approx. 100 windows
if mod(window, 2) ~= 0        % if odd
    window = window + 1;    
end
overlap = round(0.75*window); % window overlap
nfft    = 2^16;               % fft points

figure('Position', [0, 0, 0.75*screen(3), 0.75*screen(4)]);
spectrogram(x(1:4*Fs), hamming(window), overlap, nfft, Fs,...
            'yaxis','MinThreshold', -160, 'power');
set(gca,'YScale','log');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
axis([1 2 0.03 24]); view([75 55]);
set(gca,'YTick',[0.03 0.060 0.125 0.25 0.5 1 2 4 8 16 24]); 
set(gca,'YTickLabel',{30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '24k'})
title({'16-bit 1 kHz sine -105 dB peak with dither';'Spectrogram'});
c = colorbar; c.Label.String = 'Magnitude (dB)';
load('black_red_yellow.mat'); colormap(cmap);
print(gcf,'./figures/dither_16bit_1kHz_105db_spectrogram_3D','-dpng','-r150'); 

%% Normalize

if abs(ymin) >= ymax
    xnorm = x./abs(ymin);
else
    xnorm = x./ymax;
end

% export
audiowrite('./audio/dither_1kHz_105db_normalized.wav', xnorm(1:5*Fs), Fs);
