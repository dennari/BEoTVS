
w = 0.5;
q = 0.01;
r = 0.1;

% This is the transition matrix
A = [cos(w)    sin(w)/w; 
   -w*sin(w) cos(w)];

% This is the process noise covariance
Q = [0.5*q*(w-cos(w)*sin(w))/w^3 0.5*q*sin(w)^2/w^2;
   0.5*q*sin(w)^2/w^2          0.5*q*(w+cos(w)*sin(w))/w];

% SIMULATE DATA 

% This is the true initial value
x0 = [0;0.1]; 
steps = 100;
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
%figure;
%plot(T,X(1,:),'--',T,Y,'o');
methods.grid = 1;
methods.metropolis = 1;

%% find the minimum using a grid search
if methods.grid
	N1=150;
	a1 = q-q/2;
	b1 = q+3*q/2;
	X1 = linspace(a1,b1,N1);
	N2=100;
	a2 = r-r/2;
	b2 = r+r/2;
	X2 = linspace(a2,b2,N2);

	[X1,X2] = meshgrid(X1,X2);
	lh = @(q,r) energy(w,q,r,Y);
	v=zeros(N2,N1);
	ii = 0;
	pr = N1*N2/100;
	for i1=1:N1
		i1
		for i2=1:N2
			v(i2,i1)=energy(w,X1(i2,i1),X2(i2,i1),Y);
			% follow progress
			% iic = (i1-1)*N2+i2; 
			% if(iic-ii > pr)
			% 	iic/(N1*N2)
			% 	ii = iic;
			% end
		end
	end
	%[XX1,XX2] = meshgrid(X1,X2);
	mv = min(min(v));
	mmv = max(max(v));
	v_ = v;
	v(v>(mv+(mmv-mv)*0.05))=mv+(mmv-mv)*0.05;
	pcolor(X1,X2,v);
	xlabel('q');
	ylabel('r');
	q_ = X1(v==mv)
	r_ = X2(v==mv)
end

if methods.metropolis
	
end