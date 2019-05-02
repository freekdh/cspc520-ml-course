%Example 2: 
%dydx = sin(pi x)
%interval [0,1]

S = 1000;
C = rand(100,1);
D = squareform(pdist(C,'euclidean'));
%A = exp(-S*D.^2);
A = sqrt(1+S*D.^2);
B = sin(pi*C);
W = linsolve(A,B);

X = linspace(-1,2,100);
dY = arrayfun(@(z) predictdy(z,W,C,S), X);
Y = arrayfun(@(z) predicty(z,W,C,S), X);

figure(1)
plot(X,Y)
hold on 
plot(X,-cos(pi*X)/pi+1/pi+predicty(0,W,C,S))
hold off

figure(2)
plot(X,dY)
hold on 
plot(X,sin(pi*X))
hold off

function y = predicty(x,W,C,S)
y = dot(W,(sqrt(S)*(x-C).*sqrt(1+S*(x-C).^2) + asinh(sqrt(S)*(x-C)))/(2*sqrt(S)));
end

function dy = predictdy(x,W,C,S)
%dy = dot(W,exp(-S*(x-C).^2));
dy = dot(W,sqrt(1+S*(x-C).^2));
end