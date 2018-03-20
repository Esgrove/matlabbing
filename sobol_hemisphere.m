%% Sobol samples on unit hemisphere
%  Juri Lukkarila 2017

% clear everything
clc; close all; clearvars;

% screen size to position figures
screen = get(0,'screensize');

%% sobol points
N = 128;

% sobol sequence
p  = sobolset(2); 						% 2D sequence
p  = scramble(p, 'MatousekAffineOwen'); % randomize
ab = net(p,N); 							% draw #N random samples

% plot sobol values
figure('Position', [screen(3)/2-400, screen(4)/2-300, 800, 600]);
scatter(2*ab(:,1)-1, 2*ab(:,2)-1, 'filled'); grid on; 
axis([-1 1 -1 1]); axis square;

%% map to hemisphere

fov = 90;

% scaling
s = cos(0.5.*(pi - fov/180*pi));

% hemisphere for reference
[X,Y,Z] = sphere(30);

figure('Position', [screen(3)/2-600, screen(4)/2-400, 1200, 800]);
surf(X,Y,Z,'FaceAlpha',0.1, 'LineStyle','--'); grid on; hold on;
colormap([0.4 0.4 0.4]); axis([-1 1 -1 1 0 1.1]); 
for n = 1:N
	% map to disk:
    % Shirley & Chiu: "A Low Distortion Map Between Disk and Square"
    % Journal of Graphics Tools, Volume 2 Issue 3, 1997 

    a = 2 * ab(n,1) - 1; % -> [-1,1]
    b = 2 * ab(n,2) - 1;

    if (a > -b)
        if (a > b)
            r = a;
            phi = (pi / 4) * (b / a);
        else
            r = b;
            phi = (pi / 4) * (2 - (a / b));
        end

    else
        if (a < b)
            r = -a;
            phi = (pi / 4) * (4 + (b / a));    
        else
            r = -b;
            if (b ~= 0)
                phi = (pi / 4) * (6 - (a / b));
            else
                phi = 0;
            end
        end
    end
    x = s * r * cos(phi);
    y = s * r * sin(phi);
    z = sqrt(1 - x*x - y*y);
    
    % plot
    scatter3(x,y,z,'r','filled');
    
    % uncomment to visualize separately
    %pause(0.1);
end

%% View cone

[i,j,k] = cylinder([1 0],30);

h = sqrt(1 - s^2);

i = s * i;
j = s * j;
k = h * (1-k);

surf(i,j,k,'FaceAlpha',0.3, 'LineStyle','none');