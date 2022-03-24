%% SIR MODEL
%% Time interval
t0 = 0;          %initial time
tfinal = 100;    %final time
%% Initial conditions
y0=[450 450 50 50 0];  %initial sizes of susceptible1 , susceptible2, infected1, infected2 and recovered, respectively   
%% Parameters
alpha   = 0.00035; %infection rate 1
beta    = 0.01;     %infection rate 2
delta   = 0.005;     %infection rate 3
epsilon = 0.2;     %infection rate 4
rho     = 0.2;     %recovery rate
%% ODEs integration
[t,y] = ode45(@sir,[t0 tfinal],y0,[],alpha,beta,delta,epsilon, rho);
%% Plot
plot(t,y(:,1),'r',t,y(:,2),'b',t,y(:,3),'m',t,y(:,4),'g',t,y(:,5),'k','MarkerSize',3,'Linewidth',2);
xlabel('t');
xlim([0 40])
ylabel('Number of individuals');
legend('Sm(t): susceptible masked','Sn(t): susceptible unmasked','Im(t): infected masked','In(t): infected unmasked','R(t): recovered')
grid on

function Dy=sir(t,y,alpha,beta,delta,epsilon,rho)
%y(1): susceptible masked
%y(2): susceptible unmasked
%y(3): infected masked
%y(4): infected unmasked
%y(5): recovered

Dy1 = -alpha*y(1)*y(3)-delta*y(1)*y(4);
Dy2 = -beta*y(2)*y(3)-epsilon*y(2)*y(4);
Dy3 = alpha*y(1)*y(3) + delta*y(1)*y(4) - rho*y(3);
Dy4 = beta*y(2)*y(3) + epsilon*y(2)*y(4) - rho*y(4);
Dy5 = rho*y(3) + rho*y(4);

Dy=[Dy1 Dy2 Dy3,Dy4, Dy5]';
end