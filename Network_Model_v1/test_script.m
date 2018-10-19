clear all
close all
clc

X = 0:1:2500;
A1 = 2.47704;
A2 = 2.52868;
B1 = 0.00648^-1;
B2 = 0.00662^-1;
A3 = 2.49829;
B3 = .00653^-1;

Y1 = gampdf(X,A1,B1);
Y2 = gampdf(X,A2,B2);
Y3 = gampdf(X,A3,B3);

figure
plot(X,Y1,X,Y2,X,Y3)
legend('no pad','padded','one zero')

sum(X.*Y1)
sum(X.*Y2)
sum(X.*Y3)

figure

X = 50:50:1300;
c = 2/1250;

plot(X,(X-50)*-c/1250 + c)

m =  2.184e-4;

eq = @(x) x^2/m-50*x-2;
fsolve(eq,1)