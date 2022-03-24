%% SIR MODEL
%% Time interval
t0 = 0;          %initial time
tfinal = 100;    %final time
%% Initial conditions
y0=[450 450 50 50 0 900 100];  %initial sizes of susceptible1 , susceptible2, infected1, infected2, recovered,susceptible combined, and infected combined respectively   
%% Parameters
alpha   = 0.00035; %infection rate 1
beta    = 0.01;     %infection rate 2
delta   = 0.005;     %infection rate 3
epsilon = 0.2;     %infection rate 4
rho     = 0.2;     %recovery rate
%% ODEs integration
[t,y] = ode45(@sir,[t0 tfinal],y0,[],alpha,beta,delta,epsilon, rho);
%% Plot
plot(t,y(:,6),'r',t,y(:,7),'b',t,y(:,5),'MarkerSize',3,'Linewidth',2);
xlabel('t');
xlim([0 40])
ylabel('Number of individuals');
legend('S(t): susceptible', 'I(t): infected', 'R(t): recovered')
grid on

function Dy=sir(~,y,alpha,beta,delta,epsilon,rho)
%y(1): susceptible masked
%y(2): susceptible unmasked
%y(3): infected masked
%y(4): infected unmasked
%y(5): recovered
%y(6): Suceptible combined
%y(7): Infected combined

Dy1 = -alpha*y(1)*y(3)-delta*y(1)*y(4);
Dy2 = -beta*y(2)*y(3)-epsilon*y(2)*y(4);
Dy3 = alpha*y(1)*y(3) + delta*y(1)*y(4) - rho*y(3);
Dy4 = beta*y(2)*y(3) + epsilon*y(2)*y(4) - rho*y(4); 
Dy5 = rho*y(3) + rho*y(4);
Dy6 = Dy1 + Dy2;
Dy7 = Dy3 + Dy4;

Dy=[Dy1 Dy2 Dy3,Dy4, Dy5, Dy6, Dy7]';
end