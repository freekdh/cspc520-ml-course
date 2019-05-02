%Example 2: 
%dy^2dx^2 = x^2
%interval [0,1]

S = 1000;
C = rand(100,1);
D = squareform(pdist(C,'euclidean'));
%A = exp(-S*D.^2);
A = sqrt(1+S*D.^2);
B = C.*C;
W = linsolve(A,B);

X = linspace(-1,2,100);
Y = arrayfun(@(z) predicty(z,W,C,S), X);

figure(1)
plot(X,Y)
hold on 
plot(X,X.*X.*X.*X*(1/12)+predicty(0,W,C,S)+predictdy(0,W,C,S)*X)
hold off

function y = predicty(x,W,C,S)
y = dot(W,((-2+S*(x-C).^2).*sqrt(1+S*(x-C).^2)+3*sqrt(S)*(x-C).*asinh(sqrt(S)*(x-C)))/(6*S));
end

function dy = predictdy(x,W,C,S)
dy = dot(W,(sqrt(S)*(x-C).*sqrt(1+S*(x-C).^2) + asinh(sqrt(S)*(x-C)))/(2*sqrt(S)));
end

function dy2 = predictdy2(x,W,C,S)
%dy = dot(W,exp(-S*(x-C).^2));
dy2 = dot(W,sqrt(1+S*(x-C).^2));
end