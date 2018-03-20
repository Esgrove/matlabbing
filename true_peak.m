%% Sampling & True Peak Detection
%  Juri Lukkarila
%  2017

screen = get(0,'screensize'); figx = 1600; figy = 1000;
pos = [screen(3)/2-figx/2, screen(4)/2-figy/2, figx, figy];

%% read audio

[x, Fs] = audioread('get_on_it.wav');
info    = audioinfo('get_on_it.wav')

x = x(305825:397128,1);

N = length(x);              
T = 1/Fs;                      
t = 0:T:(N-1)*T;
L = max(t);

figure('Position', pos);
plot(t,x); grid on; axis([-0.02*L 1.02*L -1.05 1.05]);

%% spectrogram

window = round(N/128);      
if mod(window, 2) ~= 0      
    window = window + 1;    
end
overlap = round(0.5*window);         
nfft    = 2^nextpow2(N);             

figure('Position', pos);
spectrogram(x, hamming(window), overlap, nfft, Fs,...
            'yaxis','MinThreshold', -60, 'power');
set(gca,'YScale','log');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
axis([0 2 0.02 22]); 
set(gca,'YTick',[0.02 0.03 0.060 0.125 0.25 0.5 1 2 4 8 16 20]); 
set(gca,'YTickLabel',{20 30 60 125 250 500 '1k' '2k' '4k' '8k' '16k' '20k'})
c = colorbar; c.Label.String = 'Magnitude (dB)';
view([75 50]);

%% Oversample

r = 8;
y = resample(x, r, 1);

N2 = length(y);              
T2 = 1/(r*Fs);                      
t2 = 0:T2:(N2-1)*T2;

figure('Position', pos); hold on;
stem(t,x,'filled');
plot(t2,y); 
axis([0.002 0.006 -1.05 1.05]);

figure('Position', pos); hold on;
stem(t2,y,'Marker','none');
stem(t,x,'LineStyle','none'); 
axis([0.002 0.006 -1.05 1.05]);

%%

[x, Fs] = audioread('Bryson Tiller - Sorry Not Sorry (Esta Remix).aif');
info    = audioinfo('Bryson Tiller - Sorry Not Sorry (Esta Remix).aif')

x = x(90*Fs:150*Fs,1);

N = length(x);              
T = 1/Fs;                      
t = 0:T:(N-1)*T;
L = max(t);

figure('Position', pos);
stem(t,x,'Marker','none'); grid on; axis([-0.02*L 1.02*L -1.05 1.05]);

%% Oversample

r = 16;
y = resample(x, r, 1);

N2 = length(y);              
T2 = 1/(r*Fs);                      
t2 = 0:T2:(N2-1)*T2;

figure('Position', pos);
plot(t2,y); grid on; axis([0 60 -1.5 1.5]);

%% Plot

over = y > 1 | y < -1;

y_over  = y; y_over(~over) = NaN;

figure('Position', pos); hold on;
stem(t(34*Fs:38*Fs),x(34*Fs:38*Fs),'filled', 'MarkerSize', 3);
plot(t2(34*r*Fs:38*r*Fs),y(34*r*Fs:38*r*Fs));
plot(t2(34*r*Fs:38*r*Fs),y_over(34*r*Fs:38*r*Fs),'r');
grid on;

%%
figure('Position', [0 0 screen(3) screen(4)]);
subplot(2,2,1);
plot(t,x); grid on; axis([10 40 -1.5 1.5]); xlabel('time (s)')
ax = gca; ax.YGrid = 'on';
title('Original samples')

subplot(2,2,2);
plot(t2,y); grid on; axis([10 40 -1.5 1.5]); xlabel('time (s)')
ax = gca; ax.YGrid = 'on';
title('16x oversampling')

subplot(2,2,[3,4]);
stem(t,x,'filled', 'MarkerSize', 3); hold on;
plot(t2,y);
plot(t2,y_over,'r');
axis([35.6094 35.6136 -1.5 1.5]); xlabel('time (s)');
legend('samples', 'waveform', 'over 0 dBFS', 'Location','Best');
ax = gca; ax.YGrid = 'on';
title('True-peaks')

%%

[x, Fs] = audioread('Sorry Not Sorry -3.6db.wav');
info    = audioinfo('Sorry Not Sorry -3.6db.wav')

x = x(90*Fs:150*Fs,1);

N = length(x);              
T = 1/Fs;                      
t = 0:T:(N-1)*T;
L = max(t);

figure('Position', pos);
stem(t,x,'Marker','none'); grid on; axis([-0.02*L 1.02*L -1.05 1.05]);

% Oversample
r = 16;
y = resample(x, r, 1);

N2 = length(y);              
T2 = 1/(r*Fs);                      
t2 = 0:T2:(N2-1)*T2;

figure('Position', pos);
plot(t2,y); grid on; axis([0 60 -1.5 1.5]);

figure('Position', pos);
stem(t,x,'filled', 'MarkerSize', 3); hold on;
plot(t2,y);
axis([35.6094 35.6136 -1.5 1.5]); xlabel('time (s)');
legend('samples', 'waveform', 'Location','Best');
ax = gca; ax.YGrid = 'on';

%%

[x, Fs] = audioread('Sorry Not Sorry true-peak limit.wav');
info    = audioinfo('Sorry Not Sorry true-peak limit.wav')

x = x(90*Fs:150*Fs,1);

N = length(x);              
T = 1/Fs;                      
t = 0:T:(N-1)*T;
L = max(t);

figure('Position', pos);
stem(t,x,'Marker','none'); grid on; axis([-0.02*L 1.02*L -1.05 1.05]);

% Oversample
r = 16;
y = resample(x, r, 1);

N2 = length(y);              
T2 = 1/(r*Fs);                      
t2 = 0:T2:(N2-1)*T2;

figure('Position', pos);
plot(t2,y); grid on; axis([0 60 -1.05 1.05]);

figure('Position', pos);
stem(t,x,'filled', 'MarkerSize', 3); hold on;
plot(t2,y);
axis([35.6094 35.6136 -1.05 1.05]); xlabel('time (s)');
legend('samples', 'waveform', 'Location','Best');
ax = gca; ax.YGrid = 'on';

figure('Position', [0 0 screen(3) screen(4)]);
subplot(2,2,1);
plot(t,x); grid on; axis([10 40 -1.5 1.5]); xlabel('time (s)')
ax = gca; ax.YGrid = 'on';
title('Original samples')

subplot(2,2,2);
plot(t2,y); grid on; axis([10 40 -1.5 1.5]); xlabel('time (s)')
ax = gca; ax.YGrid = 'on';
title('16x oversampling')

subplot(2,2,[3,4]);
stem(t,x,'filled', 'MarkerSize', 3); hold on;
plot(t2,y);
axis([35.6094 35.6136 -1.5 1.5]); xlabel('time (s)');
legend('samples', 'waveform', 'over 0 dBFS', 'Location','Best');
ax = gca; ax.YGrid = 'on';
title('True-peaks')