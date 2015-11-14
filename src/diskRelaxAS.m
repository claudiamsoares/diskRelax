function [newsensors, iter] = diskRelaxAS(anchors, distances)


MYMAXITER = 100000;
epsilon = 1e-6;
nrOfAnchors = size(anchors,2);
nrOfSensors = size(distances,2)-nrOfAnchors;

MAXITER = 1e3*nrOfSensors;
dim = size(anchors,1);

%% Lipschitz constant
%% max degree
adj = distances>0;
nd = sum(adj(nrOfAnchors+1:end,nrOfAnchors+1:end),2);
dmax = max(nd);
%% max |A_i|
maxAi = max(sum(adj(1:nrOfAnchors,nrOfAnchors+1:end),1));
L = 2*dmax + maxAi;



%% function gradient
%arc-node incidence matrix
% $$$ A = triu(distances(nrOfAnchors+1:end,nrOfAnchors+1:end));
% $$$ nrOfEdges = nnz(A);
% $$$ [edgeset(:,1), edgeset(:,2)] = find(A);
% $$$ sones = ones(nrOfEdges,1);
% $$$ IXedges = (1:nrOfEdges)';
% $$$ C = kron(sparse([IXedges; IXedges],...
% $$$     edgeset,...
% $$$     [sones; -sones],...
% $$$     nrOfEdges,nrOfSensors),eye(dim));

    if isempty(anchors)
        sensors = zeros(dim,nrOfSensors);
    else
        sensors = sum(anchors,2)/nrOfAnchors*ones(1,nrOfSensors) + 0.1*randn(dim,nrOfSensors);
    end

x = sensors(:);

for node=1:nrOfSensors
    N{node} = find(distances(nrOfAnchors+node,nrOfAnchors+1:end));
end



oldx = x;


for iter=1:MAXITER
    wakeup = randi(nrOfSensors);
    update = wakeup*dim-1:(wakeup+1)*dim-2;
    for myiter=1:MYMAXITER
        w = x(update) + (myiter-2)/(myiter+1)*(x(update)-oldx(update));
        d = distances(N{wakeup}+nrOfAnchors,wakeup+nrOfAnchors);
        neighbors = reshape(x,size(sensors));
        neighbors = neighbors(:,N{wakeup});
        gradg = 0.5*calcgradh(w, neighbors, [zeros(length(N{wakeup})) d;d' 0]);
        d = distances(1:nrOfAnchors,wakeup+nrOfAnchors);
        gradh = calcgradh(w, anchors, [zeros(nrOfAnchors) d;d' 0]);
        
        gradf = gradg + gradh;
        
        oldx = x;
        x(update) = w - gradf/L;
        
        
        if norm(gradf) < epsilon
            break
        end
    end
end

newsensors = reshape(x,size(sensors));



function cost = calccostrelax(sensors, anchors, distances)
nrOfAnchors = size(anchors,2);
cost = calccostg(sensors,distances(nrOfAnchors+1:end,nrOfAnchors+1:end))...
    + calccosth(sensors, anchors, distances);


function cost = calccostg(sensors,distances)

A = triu(distances);
[edgeset(:,1), edgeset(:,2)] = find(A);

cost=0;
for m=1:size(edgeset,1)
    nodei = edgeset(m,1);
    nodej = edgeset(m,2);
    d = distances(nodei,nodej);
    cost = cost + square_pos(norm(sensors(:,nodei)-sensors(:,nodej)) - d);
end
cost = 0.5*cost;


function cost = calccosth(sensors, anchors, distances)
nrOfAnchors = size(anchors,2);
nrOfSensors = size(sensors,2);

cost = 0;
for t=1:nrOfSensors
    ancs = find(distances(1:nrOfAnchors,t+nrOfAnchors));
    for k=1:length(ancs)
        d = distances(t+nrOfAnchors,ancs(k));
        cost = cost + square_pos(norm(sensors(:,t)-anchors(:,ancs(k))) - d);
    end
end
cost = 0.5*cost;


function p = proj(x,d)

lend = length(d);
dim = length(x)/lend;
one = ones(dim,1);


z = reshape(x,[dim,lend]);
normx = sqrt(diag(z'*z));
out = one*(normx > d)';
out = out(:);
dovernormx = one*(d./normx)';
dovernormx = dovernormx(:);

p = (out.*(x.*dovernormx) ...
    + (1-out).*x);


function gradh = calcgradh(sensors, anchors, distances)
nrOfAnchors = size(anchors,2);
nrOfSensors = size(sensors,2);
dim = size(sensors,1);

gradh = zeros(dim,nrOfSensors);
for n=1:nrOfSensors
    ancs = find(distances(1:nrOfAnchors,n+nrOfAnchors));
    for m=1:length(ancs)
        a = ancs(m);
        z = sensors(:,n)-anchors(:,a);
        d = distances(n+nrOfAnchors,a);
        gradh(:,n) = gradh(:,n) + z - proj(z,d);
    end
end

gradh = gradh(:);
