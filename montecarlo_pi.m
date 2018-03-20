%% Monte Carlo integration example: Approximate pi value
%  https://en.wikipedia.org/wiki/Monte_Carlo_integration
%  Juri Lukkarila 2017

% clear everything and display more decimals
clc; close all; clearvars; format long;

% screen size to position figures
screen = get(0,'screensize');

% samples
n = 10000;

p = logspace(1,-10,n);

% draw circle
figure('Position', [screen(3)/2-400, screen(4)/2-400, 800, 800]);
phi = linspace(0,pi/2,360);
plot(sin(phi),cos(phi),'k','linewidth',1); axis square;
axis([-0.05 1.05 -0.05 1.05]); hold on; grid on;
set(gca,'Xtick', [0 1],'Ytick', [0 1]);

% run loop
m = 0; N = 0;
estimate = zeros(1,n);
error = zeros(1,n);
for i=1:n
    x = rand;
    y = rand;
    r = x^2 + y^2;
    if r < 1.0
        m = m + 1;
        plot(x,y,'b.');
    else
        plot(x,y,'r.');
    end
    value       = m/(0.25*N);
    estimate(i) = value;
    error(i)    = pi - value;
    N           = N + 1;
    
    %fprintf('pi: %.4f\n', m/(0.25*N));
    %pause(p(n));
end

fprintf('Done!\n');
fprintf('samples     = %d\n', n);
fprintf('pi estimate = %.8f\n', m/(0.25*n));
fprintf('error       = %.8f\n', abs(pi-m/(0.25*n)));

%% plot estimate and error

figure('Position', [screen(3)/2-600, screen(4)/2-300, 1200, 600]);
subplot(2,1,1);
plot(1:1:n, estimate);
axis([0 n pi-pi/6 pi+pi/6]);
set(gca,'YTick', [pi-pi/8 pi-pi/16 pi pi+pi/16 pi+pi/8]); 
set(gca,'YTickLabel',{'-\pi/8' '-\pi/16' '\pi' '+\pi/16' '+\pi/8'});  
hline = refline([0 pi]);
hline.Color = 'k';
title('Value Estimate')
grid on;

subplot(2,1,2);
plot(1:1:n, error);
axis([0 n -0.4 0.4]);
hline = refline([0 0]);
hline.Color = 'k';
title('Error');
grid on;