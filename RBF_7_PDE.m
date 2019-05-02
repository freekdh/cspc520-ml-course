%Example 2: 
%2dydx1+dydx2=sin(x1) with y(0,x2) = 0
%interval [0,1] and equidistant points

nNodes=900; %20x20
nNodesBC=20; 

S = 1000;
linsp = linspace(0,1,sqrt(nNodes));
c1 = repmat(linsp,1,sqrt(nNodes));
c2 = repelem(linsp,sqrt(nNodes));
%C = rand(nNodes,2);
C = [c1;c2].';
D = squareform(pdist(C,'euclidean'));
%A = exp(-S*D.^2);

t1=repmat(C(:,1),1,nNodes);
c1=t1.';
d1 = sqrt((t1-c1).^2);
t2=repmat(C(:,2),1,nNodes);
c2=t2.';
d2 = sqrt((t2-c2).^2);

%Interior
A1 = (1./sqrt(1+S*D.^2)).*(2*S*((2*d1)+d2));
B1 = sin(C(:,1));

%BC:
C_BC = [zeros(nNodesBC,1),linspace(0,1,nNodesBC).'];
A2 = sqrt(1+S*(C_BC.^2*ones(2,nNodes)+ones(nNodesBC,2)*(C.').^2-2*C_BC*(C.')));
B2 = zeros(nNodesBC,1);

A = [A1;A2];
B = [B1;B2];
W = linsolve(A,B);

linsp = linspace(0,1,100);
x1 = repmat(linsp,1,100);
x2 = repelem(linsp,100);
y = zeros(1,10000);

for i = 1:10000
    y(i)=predicty([x1(i),x2(i)],W,C,S);
end

figure(1)
plot3(x1,x2,y,'o')
hold on
%plot3(x1,x2,0.5*(1+2*cos(0.5*(2*x2-x1))-cos(x1)),'o')
plot3(x1,x2,0.5*(1-cos(x1)),'o')
hold off

%figure(2)
%plot3(x1,x2,y,'o')

%figure(3)
%plot3(x1,x2,0.5*(1+2*cos(0.5*(2*x2-x1))-cos(x1)),'o')

function y = predicty(x,W,C,S)
y = dot(W,sqrt(1+S*vecnorm((C-x).').^2));
end

