%% Noiseless Acoustics Audio DSP Engineer Challenge
%  Extract the fundamental frequency as a function of time
%  Juri Lukkarila / juri.lukkarila@aalto.fi / www.esgrove.fi
%  13.6.2018

clc; close all; clearvars; screen = get(0,'screensize');

%% Inport audiofile
[x, Fs] = audioread('synth_motor.wav');  
file    = audioinfo('synth_motor.wav') % print info

% normalize
x_n = x./max(abs(min(x)),max(x));

% plot signal
T = 1/Fs;                   % sample time
N = length(x);              % samples
t = 0:T:(N-1)*T;            % time vector
L = max(t);                 % length in s

% waveform
figure('Position', [screen(3)/2-1200/2, screen(4)-600, 1200, 600]);
subplot(2,1,1); plot(t,x); grid on; xlabel('time (s)');
axis([0 L -1 1]); title('Normalized audio input');

% spectrogram
window  = round(0.03*Fs);   % window
if mod(window, 2) ~= 0      % if odd
    window = window + 1;    
end
overlap = window/2;         % window overlap
nfft    = 2^15;             % fft points

subplot(2,1,2);
spectrogram(x_n, hamming(window), overlap, nfft, Fs, ...
            'yaxis','MinThreshold', -80, 'power'); axis([0 L 0 12]); 
ylabel('Frequency (kHz)'); xlabel('Time (s)'); colorbar off; 
set(gca,'XTick',0:1:9);
set(gca,'YTick',0:1:12);
axis([0 L 0 12]); 
title('Spectrogram');
print(gcf, './figures/f0_audiosignal', '-dpng', '-opengl', '-r300'); 

%% Closer look

% downsample by factor of 10: 44100 Hz -> 4410 Hz
p   = 1;
q   = 10;
x_d = resample(x_n,p,q);
Fs_d  = Fs * p / q;

% plot signal
T2 = 1/Fs_d;                % sample time
N2 = length(x_d);           % samples
t2 = 0:T2:(N2-1)*T2;        % time vector
L2 = max(t2);               % length in s

% spectrogram
window  = round(0.08*Fs_d);  % window
if mod(window, 2) ~= 0      % if odd
    window = window + 1;    
end
overlap = window/2;         % window overlap
nfft    = 2^14;             % fft points
figure('Position', [screen(3)/2-1200/2, screen(4)-600, 1200, 600]);
spectrogram(x_d, hamming(window), overlap, nfft, Fs_d, ...
            'yaxis','MinThreshold', -60, 'power'); axis([0 L2 0 2]); 
ylabel('Frequency (Hz)'); xlabel('Time (s)'); colorbar off; 
set(gca,'XTick',0:1:9);
set(gca,'YTick',0:0.2:2);
set(gca,'YTickLabel',0:200:2000);
axis([0 L2 0 2]);
title('Spectrogram');
print(gcf, './figures/f0_spectrogram_zoom', '-dpng', '-opengl', '-r300'); 

%% Linear prediction method

window_length = 0.03;             % seconds
s_full = round(window_length*Fs); % window size
s_half = round(s_full/2);         % half window size

% how many windows fits in signal with 50% overlap
windows = floor(N / s_half); windows = windows-1; index = 1;

% window function
h = hamming(s_full);

fundamentals = zeros(windows,1);
for n = 1:windows
    % windows
    x = x_n(index:index+s_full-1);
    l = length(x);
    % hamming window
    x_windowed = h .* x; 
    % scale energy to 1
    x_scaled = x_windowed./sqrt(sum(x_windowed.^2));
    % prediction order
    p = 30;
    % lpc coefficients
    b = lpc(x_scaled, p);
    % residual
    e = filter(b,1,x_scaled);
    % autocorrelation
    [ec, lags] = xcorr(e, 'coeff');
    % first positive index
    idx = find(lags >= 0, 1);
    % find largest peaks
    [peaks,p_lags] = findpeaks(ec(idx-10:end),'NPeaks',3,...
                               'SortStr','descend','MinPeakDistance',50);
    % fundamental frequency
    distance = p_lags(2)-p_lags(1);
    period = distance / Fs;
    freq = 1 / period;
    % store in list
    fundamentals(n) = freq;
    % advance index by half of window size
    index = index + s_half; 
end

figure('Position', [screen(3)/2-1200/2, screen(4)-600, 1200, 600]); hold on;
%subplot(2,1,1);
spectrogram(x_d, hamming(window), overlap, nfft, Fs_d, ...
            'yaxis','MinThreshold', -60, 'power'); axis([0 L2 0 0.8]); 
ylabel('Frequency (Hz)'); xlabel('Time (s)'); colorbar off; 
set(gca,'XTick',0:1:9);
set(gca,'YTick',0:0.1:0.8);
set(gca,'YTickLabel',0:100:800);
axis([0 L2 0 0.8]);
title('LP fundamental frequency estimation');

time = [1:1:windows].*window_length/2;
%subplot(2,1,2);
plot(time,fundamentals./1000,'LineWidth',2); hold on; grid on; 
%plot(time, smooth(fundamentals,15,'loess')); ylabel('Frequency (Hz)');
plot(time,2.*fundamentals./1000,'LineWidth',2); % second harmonic
plot(time,3.*fundamentals./1000,'LineWidth',2); % third harmonic
%axis([0 L 0 0.8]); ylabel('Frequency (Hz)'); xlabel('Time (s)');
legend('Fundamental','2nd harmonic','3rd harmonic','Location','Best');
print(gcf, './figures/f0_fundamental_lpc', '-dpng', '-opengl', '-r300'); 

%% Simple FFT-based estimation

window_length = 0.1;                % seconds
s_full = round(window_length*Fs_d); % window size
s_half = round(s_full/2);           % half window size

% how many windows fits in signal with 50% overlap
windows = floor(N2 / s_half); windows = windows-1;

% window function
h = hamming(s_full);

figure('Position', [screen(3)/2-1200/2, screen(4)-600, 1200, 600]); hold on;
spectrogram(x_d, hamming(window), overlap, nfft, Fs_d, ...
            'yaxis','MinThreshold', -60, 'power'); axis([0 L2 0 0.8]); 
ylabel('Frequency (Hz)'); xlabel('Time (s)'); colorbar off; 
set(gca,'XTick',0:1:9);
set(gca,'YTick',0:0.1:0.8);
set(gca,'YTickLabel',0:100:800);
axis([0 L2 0 0.8]);
title('FFT fundamental frequency estimation');

fundamentals = zeros(windows,1); f_index = 1;
for n = 1:windows
    % windows
    index = (n-1)*s_half + 1; 
    x = x_d(index:index+s_full-1);
    l = length(x);
    % hamming window
    x_windowed = h .* x; 
    % FFT
    nfft = 2^16;                                % fft points
    Xf = 2*abs((fft(x_windowed, nfft)/nfft));   % real part
    Xf = 20*log10(Xf./abs(max(Xf)));            % dB normalized
    f0 = Fs_d/nfft;                               % frequency resolution (Hz)
    f  = 0:f0:(nfft-1)*f0;                      % frequency vector
    index = find(f >= Fs_d/2, 1);               % index for freq Fs/2
    Xf = Xf(1:index); f = f(1:index);           % drop values over Fs/2

    % formants
    [peaks, lags] = findpeaks(Xf(find(f <= 200, 1):find(f >= 480, 1)),...
                                            'NPeaks',2,'SortStr','Descend',...
                                            'MinPeakHeight',    -30,...
                                            'MinPeakDistance',   10,...
                                            'MinPeakProminence', 20);
                                        
    if length(lags) < 2 % try again with more relaxed parameters...
        [peaks, lags] = findpeaks(Xf(find(f <= 200, 1):find(f >= 500, 1)),...
                                            'NPeaks',2,'SortStr','Descend',...
                                            'MinPeakHeight',    -60,...
                                            'MinPeakDistance',   1,...
                                            'MinPeakProminence', 1);
    end
    % simple smoothing: choose the peak that is closest to previous value
    if n > 1
        [dist, f_index] = min(abs(fundamentals(n-1)-f(lags)));
    end
    
    % fundamental frequency
    fundamentals(n) = f(lags(f_index));
    scatter(n.*window_length/2,fundamentals(n)/1000,'Filled','r');
    % spectrum
    %figure('Position',[screen(3)/2-600, screen(4)/2-300, 1200, 600]);
    %plot(f, Xf); grid on; hold on;
    %plot(f, Xf_smooth,'LineWidth', 0.25);
    %scatter(f(lags), peaks+2, 'filled', 'kv');
    % add frequency label
    %text(f(lags(1)),peaks(1)+6,strcat(num2str(round(f(lags(1)))),' Hz'));
    %text(f(lags(2)),peaks(2)+6,strcat(num2str(round(f(lags(2)))),' Hz'));
    %text(f(lags(3)),peaks(3)+6,strcat(num2str(round(f(lags(3)))),' Hz'));
    %ylabel('Magnitudi (dB)');
    %xlabel('Taajuus (Hz)');
    %legend('1/96 oktaavi', '1/9 oktaavi', 'formantit');
    %axis([0 1000 -60 0]);
    %set(gca,'XTick', 0:100:1000);
end

%print(gcf, './figures/fundamental_fft', '-dpng', '-opengl', '-r300'); 

% smoothed line
figure('Position', [screen(3)/2-1200/2, screen(4)-600, 1200, 600]); hold on;
spectrogram(x_d, hamming(window), overlap, nfft, Fs_d, ...
            'yaxis','MinThreshold', -60, 'power'); axis([0 L2 0 0.8]); 
ylabel('Frequency (Hz)'); xlabel('Time (s)'); colorbar off; 
set(gca,'XTick',0:1:9);
set(gca,'YTick',0:0.1:0.8);
set(gca,'YTickLabel',0:100:800);
axis([0 L2 0 0.8]);
title('FFT fundamental frequency estimation');

time = [1:1:windows].*window_length/2;
plot(time,smooth(fundamentals./1000),'-r','LineWidth',2); hold on; grid on; 
%axis([0 L 0 0.8]); ylabel('Frequency (Hz)'); xlabel('Time (s)');
print(gcf, './figures/f0_fundamental_fft_smooth', '-dpng', '-opengl', '-r300');