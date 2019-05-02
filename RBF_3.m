%Example 2: 
%dydx = 6*pi*(sin(2*pi*C)/(4*pi)-sin(6*pi*C)/(12*pi))
%interval [0,1]

S = 1000;
C = rand(100,1);
D = squareform(pdist(C,'euclidean'));
%A = exp(-S*D.^2);
A = sqrt(1+S*D.^2);
B = 6*pi*(sin(2*pi*C)/(4*pi)-sin(6*pi*C)/(12*pi));
W = linsolve(A,B);

X = linspace(-1,2,100);
dY = arrayfun(@(z) predictdy(z,W,C,S), X);
Y = arrayfun(@(z) predicty(z,W,C,S), X);

figure(1)
plot(X,Y,'--' ,'LineWidth',4)
hold on 
plot(X,2*(-3*cos(2*pi*X)/(8*pi)+cos(6*pi*X)/(24*pi))+2/(3*pi)+predicty(0,W,C,S),'--','LineWidth',4)
hold off

figure(2)
plot(X,dY,'--' ,'LineWidth',4)
hold on 
plot(X,6*pi*(sin(2*pi*X)/(4*pi)-sin(6*pi*X)/(12*pi)),'--' ,'LineWidth',4)
hold off

function y = predicty(x,W,C,S)
y = dot(W,(sqrt(S)*(x-C).*sqrt(1+S*(x-C).^2) + asinh(sqrt(S)*(x-C)))/(2*sqrt(S)));
end

function dy = predictdy(x,W,C,S)
%dy = dot(W,exp(-S*(x-C).^2));
dy = dot(W,sqrt(1+S*(x-C).^2));
end