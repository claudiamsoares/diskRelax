function [newsensors, varargout] = diskRelax(anchors, distances, varargin)
% DISKRELAX  Computes position estimates, given anchor
% positions and range measurements.  
%   SENSORS = DISKRELAX(A,D) where A is a matrix collecting the
%   anchor positions in the columns and D is the distance
%   matrix whose entry d_{ij} = \|x_i-x_j\| + noise if
%   measurement is available and d_{ij} = 0 otherwise; first
%   i's and j's relate to anchors.
%   SENSORS = DISKRELAX(A,D,'optname',optval) modifies the behavior of
%   the algorithm. 
%       'MAXITER': forces the algorithm to run for optval iterations;
%       'epsilon': stop when gradient norm is less than epsilon
%         [default:1e-6].
% 
%   Example: Single source localization
%   A = [0,0;0,1;1,0;1,1]'; % Anchor positions in 2D
%   X = [0.5;0.2];
%   D = dist([A, X]);
%   D(1:end-1,end) = abs(D(1:end-1,end) + 0.01*randn(size(D(1:end-1,end))));
%      
  


  [MAXITER, epsilon, stopWhenDone, savedata] = decodeVarargin(varargin);
  

  nrOfAnchors = size(anchors,2);
  nrOfSensors = size(distances,2)-nrOfAnchors;
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
  A = triu(distances(nrOfAnchors+1:end,nrOfAnchors+1:end));
  nrOfEdges = nnz(A);
  if nrOfEdges > 0
    [edgeset(:,1), edgeset(:,2)] = find(A);
  else
    edgeset = [];
  end
  sones = ones(nrOfEdges,1);
  IXedges = (1:nrOfEdges)';
  C = kron(sparse([IXedges; IXedges],...
                  edgeset,...
                  [sones; -sones],...
                  nrOfEdges,nrOfSensors),eye(dim));


  if isempty(anchors)
    sensors = zeros(dim,nrOfSensors);
  else
    sensors = sum(anchors,2)/nrOfAnchors*ones(1,nrOfSensors) + ...
              0.1*randn(dim,nrOfSensors);
  end
  
  x = sensors(:);

  
  oldx = x;
  
  if nrOfEdges > 0
    dvec = distances(sub2ind(size(distances),edgeset(:,1)+nrOfAnchors, ...
                             edgeset(:,2)+nrOfAnchors));
  else
    dvec = [];
  end

  for iter=1:MAXITER 
    w = x + (iter-2)/(iter+1)*(x-oldx);
    Cx = C*w;
    if nrOfEdges > 0
      gradg = C'* (Cx - proj(Cx,dvec));
    else
      gradg = 0;
    end
    gradh = calcgradh(reshape(w,size(sensors)), anchors, distances);
    gradf = gradg + gradh;
    
    oldx = x;
    x = w - gradf/L;
    
    if savedata.p
      [savedata] = saveIntermediateSols(iter, x, oldx, savedata);
    end
    if stopWhenDone && (norm(gradf) < epsilon)
      break
    end
  end
  
  newsensors = reshape(x,size(sensors));
  
  if nargout > 1
    varargout{1} = iter;
  end
  if nargout > 2
    varargout{2} = savedata;
  end



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
    ancs = find(distances(n+nrOfAnchors,1:nrOfAnchors)>0);
    for m=1:length(ancs)
      a = ancs(m);
      z = sensors(:,n)-anchors(:,a);
      d = distances(n+nrOfAnchors,a);
      gradh(:,n) = gradh(:,n) + z - proj(z,d);
    end
  end
  
  gradh = gradh(:);
