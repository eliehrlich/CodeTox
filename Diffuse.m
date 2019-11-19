%% Nutrient harvesting

% Simulate the movement of a number of cells that harvest from an
% underlying field of nutrients. The nutrients diffuse and are being
% replenished with linear reaction dynamics

clear all
close all

% Number of grid cells on each axis; number of cells, number of time steps
Ngrid = 100;
Ncell = 100;
Nt = 100;

% Carrying capacity of nutrients
K = 64;

% Replensihing rate of nutrients
lambda = 1;

% Diffusivity of nutrients and of cells
Dn = 1;
Dc = 1;

% Time steps
dt = 0.1;
dx = 1;
 
% Diffusion kernel (1 dimension). Warning! D/dx/dx*dt must be < 1/3
alpha = Dn/dx/dx*dt ;
if alpha > 1/3,
    error('Time step too large!')
end

kernel = alpha * [1,-2,1] + [0,1,0];

C = K*rand(Ngrid);

x = ceil(rand(Ncell,1)*Ngrid);
y = ceil(rand(Ncell,1)*Ngrid);

%% 
for t=1:Nt,
    % Diffuse columns
    for i=1:Ngrid,
        Cnew = conv(C(:,i),kernel);
        Cnew(2) = Cnew(2) + Cnew(end);
        Cnew(end-1) = Cnew(end-1) + Cnew(1);
        C(:,i) = Cnew(2:(end-1));
    end

    % Diffuse rows
    for i=1:Ngrid,
        Cnew = conv(C(i,:),kernel);
        Cnew(2) = Cnew(2) + Cnew(end);
        Cnew(end-1) = Cnew(end-1) + Cnew(1);
        C(i,:) = Cnew(2:(end-1));
    end
    
    % React (replenish) nutrients
    C = K + (C-K)*exp(-lambda*dt);
    
    % Consume
    for i=1:Ncell,
        C(x(i),y(i))=0;
    end
    
    % Move cells
    x = mod(round(x -1 + randn(Ncell,1)*2*Dc*dt),Ngrid)+1;
    y = mod(round(y -1 + randn(Ncell,1)*2*Dc*dt),Ngrid)+1;

    image(C)
    title(t)
    colorbar
    drawnow
end

%% Addendum: Solve for steady state

R = 0.1;

H = @(r) 1/2/pi/Dn*besselk(0,sqrt(lambda/Dn)*r);
[XX,YY] = meshgrid(x+rand(Ncell,1),y+rand(Ncell,1));
DD = sqrt((XX-XX').^2 + (YY-YY').^2)*dx;

for i=1:Ncell,
    DD(i,i) = R;
end

A = H(DD);

b = K*ones(Ncell,1);

q = A\b;

Cinf = K*ones(Ngrid);

for i=1:Ngrid,
    for j=1:Ngrid,
        for k=1:Ncell,
            Cinf(i,j) = Cinf(i,j) - q(k)*H( sqrt((i-x(k))^2+(j-y(k))^2)*dx);
        end
    end
end

figure
image(Cinf)
colorbar