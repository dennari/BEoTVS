%% true signal

a = -2;
b = 2;
N = 500;
d = linspace(a,b,N);
s = sin(d);%d.^3;
s_ = 3*d.^2;
plot(d,s);
%yl = ylim;
l = max(s)-min(s);
ylim([min(s)-0.1*l,max(s)+0.1*l]);
xlim([a-0.1*(b-a),b+0.1*(b-a)]);
hold on;

% measurements

std = 0.05*l;
mn = 0;
n = 150;
xi = floor(1:N/n:N);
x = d(1,xi);
y = s(1,xi)+normrnd(mn,std,1,n);
plot(x,y,'.k');
hold on;

A = ones(2); A(1,2) = x(2)-x(1);
H = [1 0];
Q = [15 2; 2 15];
R = 5*std^2;

% kalman filter

%m = [mean(y) mean(s_)]';
m = zeros(2,1);
P = eye(2);
ms = zeros(2,n);
Ps = zeros(2,2,n);
for k=1:n
    % prediction
    m_= A*m;
    P_ = A*P*A'+Q;
    % update
    S = H*P_*H'+R;
    K = P_*H'/S;
    m = m_+K*(y(k)-H*m_);
    P = P_-K*S*K';
    ms(:,k) = m;
    Ps(:,:,k) = P;
end

plot(x,ms(1,:),'-r');
plot(x,ms(2,:),'-c');
hold off;
