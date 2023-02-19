//Aim : To solve a second order linear ordinary differential equation using Euler, Modified Euler, Runge-Kutta 2nd order & 4th order methods.
//Date : 25/11/2021
clc
clear
//y=y(1)
//y'=y(2)
function dy=f(X,y)
    dy(1)=y(2) //dy(1)=y'(1)=y(2)
    dy(2)=X^3+((4*(1+X))/X)*y(2)-((2*(1+X))/X^2)*y(1)//dy(2)=y'(2)=(y')'
endfunction
X0=1
X(1)=X0
y0=[exp(2)/2;((-3*exp(2))/2)-0.5]
disp("initial condition of dependent variable : ",y0)
y1=y0
X=1:0.1:3
disp(X')
n=length(X)
h=X(2)-X(1)
//Eulers method
for i=1:n-1
    y=y0+h*f(X(i),y0)
    y0=y
    y1=[y1 y0]
end
disp("eulers : ",y1')
plot(X,y1(1,:))

//Modified Euler's method
y_M=[exp(2)/2;((-3*exp(2))/2)-0.5]
y4=y_M

for i=1:n-1
    y=y_M+h*0.5*(f(X(i),y_M)+f(X(i+1),(y_M+h*f(X(i),y_M))))
    y_M=y
    y4=[y4 y_M]
end
disp("Modified euler, y4 :",y4')
plot(X,y4(1,:),'-*m')


//R_K 2nd order
exec('C:\Users\MANIDIPA BANERJEE\Desktop\SEM III MP LAB\second order ODE\R_K  2 func for second order.sci', -1)
y_R=[exp(2)/2;((-3*exp(2))/2)-0.5]
sol_RK2=RKsecond(X0,y_R,X,f,h)
disp("RK2 solutions : ",sol_RK2')
plot(X,sol_RK2(1,:),'*r')




//R_K 4th order 
exec('C:\Users\MANIDIPA BANERJEE\Desktop\SEM III MP LAB\second order ODE\R_K 4th order for second order .sci', -1)
y_K=[exp(2)/2;((-3*exp(2))/2)-0.5]
sol_RK4=RKfourth(X0,y_K,X,f,h)
disp("RK4 solutions : ",sol_RK4')
plot(X,sol_RK4(1,:),'-^k')

//by inbuilt command
y_new=[exp(2)/2;((-3*exp(2))/2)-0.5]
sol=ode(y_new,X0,X,f)
disp("By inbuilt command : ",sol')
plot(X,sol(1,:),'-*g')
title("Plotting of Second order linear differential equation")
title color Red
title fontsize 4
xlabel("X ---->")
xlabel color magenta fontsize 4
ylabel("f(X) ----->")
ylabel color green fontsize 4
legend(["By Eulers method";"By Modified Eulers method";"By RK 2nd order";"By RK 4th order";"By inbuilt ODE"])
