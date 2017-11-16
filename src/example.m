clear all;

  %% Define sensor positions
A = [0,0;0,1;1,0;1,1]';

  %% Ground truth
X = [0.5,0.2;0.3,0.8]';

  %% Create distance matrix 
D = dist([A, X]);
D(1:4,1:4) = 0;

  %% Feed algorithms with anchors and distances
x_hat = diskRelax(A,D);
xAS_hat = diskRelaxAS(A,D);


