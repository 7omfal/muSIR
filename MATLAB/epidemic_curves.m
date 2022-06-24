%% scaled muSIR

% Population parameters
N = 1000000; %total population
s = (N-1)/N; %proportion susceptible
f = 0.70; %proportion masked

% Virus parameters
P = 3.5 %P_o = prior reproduction number
m = 0.9 %probability of mask success
p = 1 - m; %probability of mask failure
itime = 15.2 %average infected period in days
g = 1/15.2; %\gamma = recovery rate in recoveries per day
b = P*g %\beta



% Scaled time interval
t0 = 0;          %initial time
tfinal = 150;    %final time


% Integrate and plot

for i = 1:6 
    h = f + (1/10);
    f = h
    y0 = [1-f f 1/N 0 0 0];
    [t,y] = ode45(@muSIR,[t0 tfinal],y0,[],b,p,g);
    plot(t,y(:,3)+y(:,4));
    hold on
end
plot(t,y(:,3)+y(:,4));
xlabel('$\gamma t$', 'Interpreter', 'latex');
ylabel('infected fraction $i_u + i_m$', 'Interpreter','latex');
legend('$f = 0.4$', '$f = 0.5$', '$f = 0.6$', '$f= 0.7$', '$f = 0.8$','Interpreter', 'latex')
grid on

function Dy = muSIR(t,y,b,p,g)
%y(1): susceptible unmasked per total unmasked
%y(2): susceptible masked per total masked
%y(3): infected unmasked per total unmasked
%y(4): infected masked per total masked
%y(5): recovered per total population

Dy1 = -(b/g)*y(1)*(y(3) + p*y(4));
Dy2 = -(b/g)*p*y(2)*(y(3) + p*y(4));
Dy3 = (b/g)*y(1)*(y(3) + p*y(4)) - y(3);
Dy4 = (b/g)*p*y(2)*(y(3) + p*y(4)) - y(4);
Dy5 = y(3);
Dy6 = y(4);

Dy=[Dy1 Dy2 Dy3 Dy4 Dy5 Dy6]';
end
