%% scaled muSIR

% Population parameters
N = 1000; %total population
ps = 0.95; %proportion susceptible
pi = 1 - ps; %proportion infected
pm = 0.8; %proportion masked
pu = 1 - pm; %proportion unmasked

% Virus parameters
b = 0.005; %maskless infection rate
m = 0.1; %probability of mask failure
r = 0.1; %recovery rate

% Scaled time interval
t0 = 0;          %initial time
tfinal = 50;    %final time

% Initial conditions
y0 = [ps ps pi pi 0];

% Integrate and plot
[t,y] = ode45(@muSIR,[t0 tfinal],y0,[],N,pu,pm,b,m,r);
for i = 1:5    
    plot(t,y(:,i));
    hold on
end
xlabel('$\tau$', 'Interpreter', 'latex');
ylabel('scaled populations');
legend('$\hat{S}_u$', '$\hat{S}_m$', '$\hat{I}_u$', '$\hat{I}_m$', '$\hat{R}$', 'Interpreter', 'latex')
grid on

function Dy = muSIR(t,y,N,pu,pm,b,m,r)
%y(1): susceptible unmasked per total unmasked
%y(2): susceptible masked per total masked
%y(3): infected unmasked per total unmasked
%y(4): infected masked per total masked
%y(5): recovered per total population

Dy1 = -pm*y(1)*y(4) - (pu/m)*y(1)*y(3);
Dy2 = -pu*y(2)*y(3) - m*pm*y(2)*y(4);
Dy3 = pm*y(1)*y(4) + (pu/m)*y(1)*y(3) - (r/(b*N*m))*y(3);
Dy4 = pu*y(2)*y(3) + m*pm*y(2)*y(4) - (r/(b*N*m))*y(4);
Dy5 = (r/(b*N*m))*(pu*y(3) + pm*y(4));

Dy=[Dy1 Dy2 Dy3 Dy4 Dy5]';
end