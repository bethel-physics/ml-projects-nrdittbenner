clear
%%Define variables
N = 500; %number of grid points
L = 0.22;  
THb = 2.2*10^-3;
sigma = 1.2*10^-6;
kp = 5*10^-4;
km = 90;
m=3.6;
p(1) = 40;
H(1) = 1.7*10^-3;
dx = L/N; %space grid size
timesteps = input('Enter the number of time steps: ');
v = input('Enter speed of blood flow: ');
stablim = dx/v;
disp(['Stability limit = ' num2str(stablim)]) 
tau = input('Enter a time step: ');


xi = 2:N;               % Space index i
xp = 3:N+1;             % Space index i+1
xm = 1:N-1;             % Space index i-1
Time = 0;
tplot(1) = Time;
x = (0:N)*dx; 
i = 1;

%%Define Initial Conditions
P = zeros(N+1,1);
H = zeros(N+1,1);
P(1:N/2) = 44;
P(N/2:N) = 102;
H(1:N/2) = 0;
H(N/2:N) = 0.5;
%%Fill matrix
method = input('Choose a method (1 = original, 2 = FTCS, 3 = Lax Wenroff,4 = Upwind): ')
option = input('With constants (1 for yes, 2 for no)? ');
coef = (v*tau)/(2*dx);
coef2 = (v^2*tau^2)/(2*dx^2);
while (i <= timesteps)
    Time = Time +tau;
    tplot(i+1) = Time;
    Pplot(i+1,:) = P;
    HPlot(i+1,:) = H;
    if method == 1
        if option == 2
            P(xi) =tau*(-v*((P(xp)-P(xm))/(2*dx)))+P(xi);
            H(xi) = tau*(-v*((H(xp)-H(xm))/(2*dx)))+H(xi);
        else
            P(xi) =tau*(-v*((P(xp)-P(xm))/(2*dx))+(m*((km.*H(xi)-kp.*(THb-H(xi))).*P(xi))/sigma))+P(xi);
            H(xi) = tau*(-v*((H(xp)-H(xm))/(2*dx))+kp.*(THb-H(xi)).*P(xi)-km.*H(xi))+H(xi);
        end
    elseif method == 2
    %FTCS
        if option == 2
            P(xi) = P(xi)-((v*tau)/(2*dx))*(P(xp)-P(xm));
            H(xi) = H(xi)-((v*tau)/(2*dx))*(H(xp)-H(xm));
        else
            P(xi) = P(xi)-((tau*v/(2*dx))*(P(xp)-P(xm))+(m*((km.*H(xi)-kp.*(THb-H(xi))).*P(xi))/sigma));
            H(xi) = H(xi)-((v*tau/(2*dx))*(H(xp)-H(xm))+kp.*(THb-H(xi)).*P(xi)-km.*H(xi)); 
        end
    elseif method == 3
    %Lax-Wenroff
        if option == 1
            P(xi) = P(xi)-coef*(P(xp)-P(xm))+coef2*(P(xp)-2*P(xi)+P(xm))+(m*((km.*H(xi)-kp.*(THb-H(xi))).*P(xi))/sigma);
            H(xi) = H(xi)-coef*(H(xp)-H(xm))+coef2*(H(xp)-2*H(xi)+H(xm))+kp.*(THb-H(xi)).*P(xi)-km.*H(xi);
        else
            P(xi) = P(xi)-coef*(P(xp)-P(xm))+coef2*(P(xp)-2*P(xi)+P(xm));
            H(xi) = H(xi)-coef*(H(xp)-H(xm))+coef2*(H(xp)-2*H(xi)+H(xm)); 
        end
    elseif method == 4
        if option == 2
            P(xi) = P(xi)-(tau/dx)*(v*(P(xi)-P(xm)));
            H(xi) = H(xi)-(tau/dx)*(v*(H(xi)-H(xm)));
        else
            P(xi) = P(xi)-(tau/dx)*(v*(P(xi)-P(xm)))+(m*((km.*H(xi)-kp.*(THb-H(xi))).*P(xi))/sigma);
            H(xi) = H(xi)-(tau/dx)*(v*(H(xi)-H(xm)))+kp.*(THb-H(xi)).*P(xi)-km.*H(xi);
        end
    end
    i = i+1;
end
surf(x,tplot,Pplot,'EdgeColor','none');xlabel('space');ylabel('time');zlabel('pressure');