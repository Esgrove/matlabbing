%% Monte Carlo integration example: Approximate pi value
%  https://en.wikipedia.org/wiki/Monte_Carlo_integration
%  Juri Lukkarila 2017

% clear everything and display more decimals
clc; close all; clearvars; format long;

% screen size to position figures
screen = get(0,'screensize');

%% Uniform sampling

% samples
N = 1200;

X = rand(1,N);
Y = rand(1,N);
R = X.^2 + Y.^2;

idx_in  = find(R  < 1.0);
idx_out = find(R >= 1.0);

estimate = length(idx_in)/(0.25*N);
error    = pi - estimate;
fprintf('Uniform sampling\n');
fprintf('...samples     = %d\n', N);
fprintf('...pi estimate = %.8f\n', estimate);
fprintf('...error       = %.8f\n', abs(error));

% plot
figure('Position', [screen(3)/2-400, screen(4)/2-300, 800, 600]);
phi = linspace(0,pi/2,360); % draw circle
plot(sin(phi),cos(phi),'k','linewidth',1); axis square;
axis([-0.05 1.05 -0.05 1.05]); hold on; grid on;
set(gca,'Xtick', [0 1], 'Ytick', [0 1]);
title('Pi value approximation using random sampling');
xlabel('x'); ylabel('y');
text(0.4,-0.025, join(['\pi estimate = ' sprintf('%.5f',estimate)]));
% plot samples
plot(X(idx_in ),Y(idx_in ),'.b'); % blue dots
plot(X(idx_out),Y(idx_out),'.r'); %  red dots
print(gcf,'./figures/montecarlo_pi_uniform', '-dpng', '-r300');

%% Sobol sampling

% sobol sequence
p  = sobolset(2); 						% 2D sequence
p  = scramble(p, 'MatousekAffineOwen'); % randomize
ab = net(p,N); 							% draw #N random samples

X = ab(:,1);
Y = ab(:,2);

R = X.^2 + Y.^2;

idx_in  = find(R  < 1.0);
idx_out = find(R >= 1.0);

estimate = length(idx_in)/(0.25*N);
error    = pi - estimate;
fprintf('\nSobol sampling\n');
fprintf('...samples     = %d\n', N);
fprintf('...pi estimate = %.8f\n', estimate);
fprintf('...error       = %.8f\n', abs(error));

% plot sobol values
figure('Position', [screen(3)/2-400, screen(4)/2-300, 800, 600]);
plot(sin(phi),cos(phi),'k','linewidth',1); axis square; hold on; grid on;
axis([-0.05 1.05 -0.05 1.05]); set(gca,'Xtick', [0 1], 'Ytick', [0 1]);
title('Pi value approximation using random sampling, Sobol samples'); 
xlabel('x'); ylabel('y');
text(0.4,-0.025, join(['\pi estimate = ' sprintf('%.5f',estimate)]));
plot(X(idx_in ),Y(idx_in ),'.b'); % blue dots
plot(X(idx_out),Y(idx_out),'.r'); %  red dots
print(gcf,'./figures/montecarlo_pi_sobol', '-dpng', '-r300');
fprintf('Done!\n');

%% Animation

p = logspace(1,-10,N); % decreasing pause time vector

% draw circle
figure('Position', [screen(3)/2-400, screen(4)/2-300, 800, 600]);
phi = linspace(0,pi/2,360);
plot(sin(phi),cos(phi),'k','linewidth',1); axis square;
axis([-0.05 1.05 -0.05 1.05]); hold on; grid on;
set(gca,'Xtick', [0 1], 'Ytick', [0 1]);
title('Pi value approximation using random sampling');
xlabel('x'); ylabel('y');

v = VideoWriter('./figures/montecarlo_pi_animation.mp4','MPEG-4'); 
v.FrameRate = 60;
open(v);

frameskip = 5; % how many "frames" to skip between animation

% run loop
m = 0; value = 0.0;
t = text(0.4,-0.025, join(['\pi estimate = ' num2str(round(value,6))]));
estimate = zeros(1,N);
error    = zeros(1,N);
frame    = struct('cdata', cell(1,floor(N/frameskip)),...
               'colormap', cell(1,floor(N/frameskip)));
fprintf('\nAnimation...\n');
for n = 1:N
    x = X(n);
    y = Y(n);
    r = x^2 + y^2;
    if r < 1.0
        m = m + 1;
        plot(x,y,'b.');
    else
        plot(x,y,'r.');
    end
    value       = m/(0.25*n);
    estimate(n) = value;
    error(n)    = pi - value;
    
    %fprintf('pi: %.4f\n', value));
    if mod(n, frameskip) == 0     % pause every "frameskip" frame
        pause(p(N));                % use pause to animate
        set(t,'String',join(['\pi estimate = ' sprintf('%.5f',value)]));
        frame(n/frameskip) = getframe(gcf);   % store image as videoframe
    end
end

fprintf('Done!\n');
fprintf('...samples     = %d\n', N);
fprintf('...pi estimate = %.8f\n', m/(0.25*N));
fprintf('...error       = %.8f\n', abs(pi-m/(0.25*N)));

% export video
writeVideo(v,frame);
close(v);
fprintf('Video export done!\n');

%% plot estimate and error

figure('Position', [screen(3)/2-600, screen(4)/2-300, 1200, 600]);
subplot(2,1,1);
plot(1:1:N, estimate);
axis([0 N pi-pi/6 pi+pi/6]);
set(gca,'YTick', [pi-pi/8 pi-pi/16 pi pi+pi/16 pi+pi/8]); 
set(gca,'YTickLabel',{'-\pi/8' '-\pi/16' '\pi' '+\pi/16' '+\pi/8'});  
hline = refline([0 pi]);
hline.Color = 'k';
title('Value Estimate'); xlabel('sample (n)'); grid on;

subplot(2,1,2);
plot(1:1:N, error);
axis([0 N -0.4 0.4]);
hline = refline([0 0]);
hline.Color = 'k';
title('Error'); xlabel('sample (n)'); grid on;
print(gcf,'./figures/montecarlo_pi_error', '-dpng', '-r300');