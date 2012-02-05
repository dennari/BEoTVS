  randn('state',123);
  steps = 100;

  w = 0.5;
  q = 0.01;
  r = 0.1;

  % This is the transition matrix
  A = [cos(w)    sin(w)/w; 
       -w*sin(w) cos(w)];

  % This is the process noise covariance
  Q = [0.5*q*(w-cos(w)*sin(w))/w^3 0.5*q*sin(w)^2/w^2;
       0.5*q*sin(w)^2/w^2          0.5*q*(w+cos(w)*sin(w))/w];

  % This is the true initial value
  x0 = [0;0.1]; 

  X = zeros(2,steps);  % The true signal
  Y = zeros(1,steps);  % Measurements
  T = 1:steps;         % Time
  x = x0;
  for k=1:steps
    x = mvnrnd(A*x,Q)';
    y = mvnrnd(x(1),r);
    X(:,k) = x;
    Y(:,k) = y;
  end

 figure;
 subplot(1,1,1);
 plot(T,X(1,:),'--',T,Y,'o');
 legend('True signal','Measurements');