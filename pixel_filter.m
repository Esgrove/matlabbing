%% Filter shapes
%  Juri Lukkarila
%  2017

clc; close all; clearvars; screen = get(0,'screensize');

%% Mitchell-Netravali 1D

N = 100;                % points
x = linspace(-2, 2, N); % vector
y = zeros(1,N);         % pre-allocate result
for n = 1:N
    y(n) = MitchellNetravali(x(n));
end

figure('Position', [screen(3)/2-400, screen(4)/2-300, 800, 600]);
plot(x,y); grid on; axis([-2 2 -0.1 1.1]); set(gca,'YTick',0:0.1:1); 
title('Mitchell-Netravali 1D');

%% Mitchell-Netravali 2D

N = 64;

t = linspace(-2, 2, N);

[X,Y] = meshgrid(t);

R = sqrt(X.^2 + Y.^2) + eps;

Z = zeros(N,N);
for n = 1:N
    for m = 1:N
        Z(n,m) = MitchellNetravali(R(n,m));
    end
end

% grid mesh
figure('Position', [screen(3)/2-800, screen(4)/2-300, 800, 600]);
surf(X,Y,Z); grid on; axis([-2 2 -2 2 -0.1 1.1]); set(gca,'ZTick',0:0.1:1);
title('Mitchell-Netravali 2D');

%% Mitchell-Netravali 2D shaded

N = 256;

t = linspace(-2, 2, N);

[X,Y] = meshgrid(t);

R = sqrt(X.^2 + Y.^2) + eps;

Z = zeros(N,N);
for n = 1:N
    for m = 1:N
        Z(n,m) = MitchellNetravali(R(n,m));
    end
end

figure('Position', [screen(3)/2, screen(4)/2-300, 800, 600]);
surf(X,Y,Z,'FaceColor','red','EdgeColor','none'); grid on; 
axis([-2 2 -2 2 -0.1 1.1]); set(gca,'ZTick', 0:0.1:1);
title('Mitchell-Netravali 2D');
camlight left; lighting phong; % shading

%% Mitchell-Netravali x*y

N = 96;

t = linspace(-2, 2, N);

[X,Y] = meshgrid(t);

Z = zeros(N,N);
for n = 1:N
    for m = 1:N
        Z(n,m) = MitchellNetravali(X(n,m)) * MitchellNetravali(Y(n,m));
    end
end

% grid mesh
figure('Position', [screen(3)/2-800, screen(4)/2-300, 800, 600]);
surf(X,Y,Z); grid on; axis([-2 2 -2 2 -0.1 1.1]); set(gca,'ZTick',0:0.1:1); 
title('Mitchell-Netravali x*y');

%% Gauss 1D

x = linspace(-2, 2, 1000);

sigma = 0.4;

y = 1 / (sqrt(2*pi) * sigma) .* exp(-x.^2 / (2 * sigma^2));

figure('Position', [screen(3)/2-400, screen(4)/2-300, 800, 600]);
plot(x,y); grid on; axis([-2 2 -0.1 1.1]); set(gca,'YTick',0:0.1:1); 
title('Gauss 1D'); hold on;

s = 0.5:0.1:0.8;

for n = 1:length(s)
    sigma = s(n);
    y = 1 / (sqrt(2*pi) * sigma) .* exp(-x.^2 / (2 * sigma^2));
    plot(x,y);
end

legend('0.4','0.5','0.6','0.7','0.8');

%% Gauss 2D

N = 64;

t = linspace(-2, 2, N);

sigma = 0.5;

[X,Y] = meshgrid(t);

Z = 1 / (sqrt(2*pi) * sigma) .* exp(-1* (X.^2 + Y.^2) / (2 * sigma^2));

figure('Position', [screen(3)/2-800, screen(4)/2-300, 800, 600]);
surf(X,Y,Z); grid on; axis([-2 2 -2 2 -0.1 1.1]); set(gca,'ZTick',0:0.1:1);
title('Gauss 2D');

%% Gauss 2D shaded

N = 512;

t = linspace(-2, 2, N);

[X,Y] = meshgrid(t);

Z = 1 / (sqrt(2*pi) * sigma) .* exp(-1* (X.^2 + Y.^2) / (2 * sigma^2));

figure('Position', [screen(3)/2, screen(4)/2-300, 800, 600]);
surf(X,Y,Z,'FaceColor','red','EdgeColor','none'); grid on;
axis([-2 2 -2 2 -0.1 1.1]); set(gca,'ZTick',0:0.1:1);
lighting phong; camlight left; title('Gauss 2D');

%% Gauss x*y

N = 64;

t = linspace(-2, 2, N);

[X,Y] = meshgrid(t);

Z = zeros(N,N);
for n = 1:N
    for m = 1:N
        Z(n,m) = Gauss(X(n,m)) * Gauss(Y(n,m));
    end
end

% grid mesh
figure('Position', [screen(3)/2-800, screen(4)/2-300, 800, 600]);
surf(X,Y,Z); grid on; axis([-2 2 -2 2 -0.1 1.1]); set(gca,'ZTick',0:0.1:1); 
title('Gauss x*y');

%% Filter functions

function [ y ] = MitchellNetravali( x )
% MitchellNetravali
%
	B = 1 / 3;
    C = 1 / 3;
	ax = abs(x);
	if (ax < 1)
		y = ((12- 9 * B - 6 * C) * ax * ax * ax + (-18 + 12 * B + 6 * C) * ax * ax + (6 - 2 * B)) / 6;
    elseif (ax < 2)
		y = ((-B - 6 * C) * ax * ax * ax + (6 * B + 30 * C) * ax * ax + (-12 * B - 48 * C) * ax + (8 * B + 24 * C)) / 6;
    else
		y = 0;
    end
end

function [ y ] = Gauss( x )
% Gaussian filter
%
sigma = 0.5;

y = 1 / (sqrt(2*pi) * sigma) .* exp(-x.^2 / (2 * sigma^2));

end