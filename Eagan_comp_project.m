%QUESTION 1 AND 2
d = 47; %diameter in mm
d = d/1000; %Diameter of cylinder in m
l = 100; %length in mm
l = l/1000; %length of cylinder in m
den = 2700; %Denisty of aluminum in kg/m^3
c = 903; %Specific heat of aluminum in J/(kg*C)
k = 241; %Thermal conductivity of aluminum = 241 W/(m*C)
A = 2*pi*(d/2)*(l+(d/2)); %area of cylinder
m = pi*(d/2)^2*l*den; %mass of cylinder

t = [0 5 10 15 20 25 30 35 40]; %time in seconds
Ta = [23.2 22.5 17.6 16.4 15.2 14.4 13.8 13.2 12.7]; %Temperature of the aluminum cylinder in C
Tw = 3.8; %Temperature of the ice water bath in C

To = Ta(1,1); %intial temperature of cylinder

%QUESTION 3 AND 4
% Newton-Raphson 3rd data point
syms h
datapoint = 3;
temp = Ta(1,datapoint); %3rd data point value
h_prev3 = 500;
t_t = Tw + (To - Tw)*exp(-(h*A*t(1,datapoint))/(m*c)) -temp; %equation 4
t_td = diff(t_t,h); %derivative of equation 4
tol = 0.02;
error = 1;
while error > tol
    f = double(subs(t_t,h,h_prev3));
    fd = double(subs(t_td,h,h_prev3));
    h_new3 = h_prev3 - (f/fd); %newton raphson equation
    error = abs((h_new3 - h_prev3)/h_new3); 
    h_prev3 = h_new3;
end

%newton rasphon datapoint 6 
syms h
datapoint = 6;
temp = Ta(1,datapoint); %3rd data point value
h_prev6 = 500;
t_t = Tw + (To - Tw)*exp(-(h*A*t(1,datapoint))/(m*c)) -temp; %equation 4
t_td = diff(t_t,h); %derivative of equation 4
tol = 0.02;
error = 1;
while error > tol
    f = double(subs(t_t,h,h_prev6));
    fd = double(subs(t_td,h,h_prev6));
    h_new6 = h_prev6 - (f/fd); %newton raphson equation
    error = abs((h_new6 - h_prev6)/h_new6); 
    h_prev6 = h_new6;
end

%Newton Rasphon data point 9
syms h
datapoint = 9;
temp = Ta(1,datapoint); %3rd data point value
h_prev9 = 500;
t_t = Tw + (To - Tw)*exp(-(h*A*t(1,datapoint))/(m*c)) -temp; %equation 4
t_td = diff(t_t,h); %derivative of equation 4
tol = 0.02;
error = 1;
while error > tol
    f = double(subs(t_t,h,h_prev9));
    fd = double(subs(t_td,h,h_prev9));
    h_new9 = h_prev9 - (f/fd); %newton raphson equation
    error = abs((h_new9 - h_prev9)/h_new9); 
    h_prev9 = h_new9;
end

%QUESTION 5
%Differentiation foward difference
deltat = 5; %seconds
n = length(Ta); %number of element in matrix
fprime = zeros(1,n-1); %store answers in this matrix
for i = 1:n-1
    fprime(1,i) = (Ta(1,i+1) - Ta(1,i))/deltat; %forward difference equation
end

%QUESTION 6
%multiplied by mc
ftime = zeros(1,n-1);
fTa = zeros (1,n-1);
for j = 1:n-1
    ftime(1,j) = t(1,j);
    fTa(1,j) = Ta(1,j);
end

vector = m*c*fprime;
figure
plot(ftime, vector) %plot of time vs differentiation multiplied by mc
hold on
xlabel('Time(s)')
ylabel('M*C*derivative')


h_avg = (h_prev3 + h_prev6 + h_prev9)/3; %averge h vvalue
heat = -h_avg*A*(fTa - Tw);
plot(ftime, heat)
legend('mc(dT/dt)', '-hA(T-Tw)')

%Question 7
%regression
p = polyfit(t,Ta,2);
a2 = p(1);
a1 = p(2);
a0 = p(3);

ta_fit = a0 + a1*t + a2*t.^2; %polynomial equation
figure
plot(t, ta_fit, '-')
hold on
plot (t, Ta,'s' )

xlabel('Time(s)')
ylabel('Temperature(C)')
legend('Second-order polynomial', 'Experimental data')
title('Second-order polynomial regression')

% QUESTION 8
% trapazoidal integration
value = length(t);
temp_diff = Ta -Tw; %temperature differene vector
 
a = t(1); %left point
b = t(value); %right point
segm = value - 1; %number of segments 
htrap = (b-a)/segm;
total = 0;
for tr = 1:segm
    total = total + (htrap/2)*(temp_diff(tr)+temp_diff(tr+1));
end
resultq8 = h_avg*A*total;
disp(resultq8)

%Question 9 
%Differential first order 
dt = 0.5; %detla t in seconds
t_values = 0:dt:40;
u_euler = zeros(1,length(t_values)); %u answer vector set up
u_heun = zeros(1,length(t_values));
u_euler(1)= To; %intial u euler value
u_heun(1) = To; %intial u heun value
syms u
g= -(h_avg*A/(m*c))*(u-Tw);
%euler method
 
for e = 1:length(t_values)-1
    g_value = double(subs(g,u,u_euler(e)));
    u_euler(e+1) = u_euler(e) + g_value*dt;
end
% Heun method 

for hm = 1:length(t_values)-1
    k1 = double(subs(g,u,u_heun(hm)));
    k2 = double(subs(g,u,u_heun(hm)+k1*dt));
    u_heun(hm+1) = u_heun(hm) + ((1/2)*k1 + (1/2)*k2)*dt;
end

u_exact = Tw + (To-Tw)*exp(-(h_avg*A/(m*c))*t_values); %exact solution
figure %plot
plot(t_values,u_euler, 'o') 
hold on
plot(t_values,u_heun, 's')
plot(t_values,u_exact, '-')
xlabel('Time(s)')
ylabel('Temperature (C)')
legend('Euler method','Heun method','Exact solution')
 
%Answer displays
fprintf('\nQuestion 2\n')
fprintf('Cylinder diameter = %.4f m\n',d)
fprintf('Cylinder length = %.4f m\n',l)
fprintf('\nQuestions 3 and 4\n')
fprintf('h using data point 3 = %.4f\n', h_prev3)
fprintf('h using data point 6 = %.4f\n', h_prev6)
fprintf('h using data point 9 = %.4f\n', h_prev9)
fprintf('Average h = %.4f\n', h_avg)
fprintf('\nQuestion 5\n')
fprintf('fprime = %.4f\n',fprime)
fprintf('\nQuestion 7\n')
fprintf('a0 = %.4f\n',a0)
fprintf('a1 = %.4f\n',a1)
fprintf('a2 = %.4f\n',a2)
fprintf('\nQuestion 8\n')
fprintf('Released energy = %.4f\n',resultq8)
fprintf('\nQuestion 9\n')
fprintf('Euler u = %.4f\n', u_euler(end))
fprintf('Heun u = %.4f\n',u_heun(end))
fprintf('Exact temperature = %.4f\n',u_exact(end))
