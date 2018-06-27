%% Audio blur
%  Juri Lukkarila
%  2018

% figure positioning
screen = get(0,'screensize'); figx = 1200; figy = 600; clc; close all;

%% Audio sample

[x, Fs] = audioread("This Is A Journey Into Sound.aif");
audio1  = audioinfo("This Is A Journey Into Sound.aif") % print info

%soundsc(x,Fs);

% use only one channel and normalize
x_n = 0.99.*(x(:,1)./max(abs(min(x(:,1))),max(x(:,1)))); x_n = x_n';

audiowrite('./audio/audioblur_in1.wav',x_n,Fs);

%% Blur

L = [2 1 1.5 0.5 0.25]; % window lengths (s) 

% batch process with different window lengths (in parallel)
for l = 1:length(L)
    n = round(Fs * L(l)); % samples
    
    % blur kernel
    %env = hamming(n);
    env = gausswin(n, 3);
    env = env - min(env); % zero and normalize to 1
    env = env./max(env);

    kernel = env .* rand(n,1); kernel = kernel';

    %figure('Position',[screen(3)/2-figx/2, screen(4)/2-figy/2, figx, figy]); 
    %stem(kernel,'Marker','None');
    %title(strcat('Blur kernel, length: ',num2str(L(l)),'s')); 
    %xlabel('sample (n)');
    %axis([-n/20 n+n/20 0 1]);

    % convolve with blur kernel
    blur = conv(x_n, kernel); 
    blur = 0.99.*(blur./max(abs(min(blur)),max(blur))); % normalize
    
    % export audio
    audiowrite(strcat('./audio/audioblur_in1_blurred_',...
        num2str(L(l)),'.wav'),blur,Fs);

    % Plot
    x = [zeros(1,round(n/2)) x_n]; % pad input by n/2

    T  = 1/Fs;

    N1 = length(x);  t1 = 0:T:(N1-1)*T; 
    N2 = length(blur); t2 = 0:T:(N2-1)*T;
    
    lx = max(max(t1),max(t2));
    
    figure('Position', [screen(3)/2-figx/2, screen(4)/2-figy/2, figx, figy]);
    subplot(2,1,1);
    plot(t1,x); grid on; axis([0 lx -1 1]); xlabel('time (s)');
    title('Original audio');
    subplot(2,1,2);
    plot(t2,blur); grid on; axis([0 lx -1 1]); xlabel('time (s)');
    text = [' blur length ', num2str(L(l)), ' s'];
    title(strcat('Blurred audio,',text));
    numstr = strrep(num2str(L(l)),'.',''); % remove dots from number
    print(gcf,strcat('./figures/audioblur_in1_',numstr), '-dpng', '-r300');
    close all;
    fprintf('L = %g done!\n',L(l)); 
end

%% Audio sample 2

[x, Fs] = audioread("cafe.wav");
audio2 = audioinfo("cafe.wav")

%soundsc(x,Fs);

% use only one channel and normalize
x_n = 0.99.*(x(:,1)./max(abs(min(x(:,1))),max(x(:,1)))); x_n = x_n';

audiowrite('./audio/audioblur_in2.wav',x_n,Fs);

%% Blur 2

L = [2 1 1.5 0.5 0.25]; % window lengths (s) 

% batch process with different window lengths (in parallel)
for l = 1:length(L)
    n = round(Fs * L(l)); % samples
    
    % blur kernel
    %env = hamming(n);
    env = gausswin(n, 3);
    env = env - min(env); % zero and normalize to 1
    env = env./max(env);

    kernel = env .* rand(n,1); kernel = kernel';

    %figure('Position',[screen(3)/2-figx/2, screen(4)/2-figy/2, figx, figy]); 
    %stem(kernel,'Marker','None');
    %title(strcat('Blur kernel, length: ',num2str(L(l)),'s')); 
    %xlabel('sample (n)');
    %axis([-n/20 n+n/20 0 1]);

    % convolve with blur kernel
    blur = conv(x_n, kernel); 
    blur = 0.99.*(blur./max(abs(min(blur)),max(blur))); % normalize
    
    % export audio
    audiowrite(strcat('./audio/audioblur_in2_blurred_',...
                num2str(L(l)),'.wav'),blur,Fs);

    % Plot
    x = [zeros(1,round(n/2)) x_n]; % pad input by n/2

    T  = 1/Fs;

    N1 = length(x);  t1 = 0:T:(N1-1)*T; 
    N2 = length(blur); t2 = 0:T:(N2-1)*T;
    
    lx = max(max(t1),max(t2));
    
    figure('Position', [screen(3)/2-figx/2, screen(4)/2-figy/2, figx, figy]);
    subplot(2,1,1);
    plot(t1,x); grid on; axis([0 lx -1 1]); xlabel('time (s)');
    title('Original audio');
    subplot(2,1,2);
    plot(t2,blur); grid on; axis([0 lx -1 1]); xlabel('time (s)');
    text = [' blur length ', num2str(L(l)), ' s'];
    title(strcat('Blurred audio,',text));
    numstr = strrep(num2str(L(l)),'.',''); % remove dots from number
    print(gcf,strcat('./figures/audioblur_in2_',numstr), '-dpng', '-r300');
    close all;
    fprintf('L = %g done!\n',L(l)); 
end