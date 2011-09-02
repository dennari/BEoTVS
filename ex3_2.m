% for outputting figures correctly
latexTextW = 341; % in points
figAspectRatio = 16/10; % w/h
figW = 1.3*latexTextW;
figH = figW/figAspectRatio;

%% Exercise round 3, exercise 2

%true signal
N=100;
q=3;

s=cumsum([0 normrnd(0,sqrt(q),1,N-1)]);
figure(1);
clf;
true_m = plot(1:N,s);
xlabel('Time');
ylabel('Mean');
hold on;

% measurements
A = 1;
H = 1;
r = 5;
mn = 0;
n = N;
x = floor(1:N/n:N);
y = s(1,x)+normrnd(mn,sqrt(r),1,n);
measurements = plot(x,y,'.k');

% Kalman filter
m = 0;
P = 1;
m2 = 0;
P2 = 1;
ms = zeros(1,n);ms(1)=m;
ms2 = zeros(1,n);
Ps2 = zeros(1,n);
Ps = zeros(1,n);Ps(1)=P;
Ks = Ps;
for k=2:n
    % prediction
    m_= m;
    P_ = P+q;
    % update
    K = P_/(P_+r);
    m = m_+K*(y(k)-m_);
    P = P_-(P_^2/(P_+r));
    ms(k) = m;
    Ps(k) = P;
    Ks(k) = K;
end
figure(1);
kalman_m = plot(x,ms,'-r');
hold on;
mso = ms;
Pso = Ps;
rmse(x,ms)
figure(2);
kalman_P = plot(x,Ps,'-r');
xlabel('Time');
ylabel('Variance');
hold on;

% discretization
a=min(y);
b=max(y);
N = 500;
t = linspace(a,b,N);
[TX,TY] = meshgrid(t);
p_dyn = normpdf(TX,TY,sqrt(q));
m = 0;
P = 1;
p_ = normpdf(t,m,P);
ms = zeros(1,n);ms(1)=m;
Ps = zeros(1,n);Ps(1)=P;
distr = zeros(N,n);
for k=2:n
    p = sum(p_dyn.*repmat(p_',1,N));
    p = normpdf(y(k),t,sqrt(r)).*p;
    p = p/sum(p);
    p_ = p;
    distr(:,k) = p';
    m = t*p';
    ms(k) = m; % mean
    Ps(k) = (t-m).^2*p';
end;

%plot
figure(1);
discr_m = plot(x,ms,'-g','LineWidth',2.0);
legend([true_m measurements kalman_m discr_m],'True','Measurements','Kalman','Discretization');
exportplot('ex_3_2_means.pdf',figW,figH);
figure(2);
discr_P = plot(x,Ps,'-g');
legend([kalman_P discr_P],'Kalman','Discretization');
set(gcf,'Color','none','Units','points','Position',[400 400 figW figH]);
exportplot('ex_3_2_variances.pdf',figW,figH);



%figure(3);
%pcolor(1:n,t,distr);
%shading interp;
%ylim([a b]);


%rmse(x,ms)/rmse(x,mso)


%% stationary Kalman filter (use the last value of K (the Kalman gain) from
% the normal filter)

m = 0;
P = 1;
ms = zeros(1,n);
Ps = zeros(1,n);
for k=1:n
    m = (1-K)*m+K*y(k);
    P = P-K^2/(P+r);
    ms(k) = m;
    Ps(k) = P;
end
kalmansolution = plot(x,ms,'-k');
ms_st = ms;
Ps_st = Ps;
rmse(x,ms)




%% RTS smoother

mss = zeros(1,n);
Pss = zeros(1,n);
Gs = mss;
m = mso(end);
P = Pso(end);
mss(end) = m;
Pss(end) = P;

for k=n-1:-1:1
    g = Pso(k)/(Pso(k)+q);
    m = mso(k)+g*(m-mso(k));
    P = Pso(k)+g^2*(P-Pso(k)-q);
    Gs(k) = g;
    mss(k) = m;
    Pss(k) = P;
end

figure(1);
rtssolution = plot(x,mss,'-g','LineWidth',2);
rmse(s,mss)

hold on;
figure(2);
plot(x,Pss,'-g','LineWidth',2);
hold on;
%yl = ylim;
%ylim([0,yl(2)]);
%figure(3);
%plot(x,mss,'-g','LineWidth',2);


if 0
% RTS smoother -- discrete

mss = zeros(1,n);
Pss = zeros(1,n);
m = mso(end);
P = Pso(end);
p_ = normpdf(t,m,P);
size(p_)
distr = zeros(N,n);
for k=n-1:-1:1
    % dynamical distribution times the "previous" smoothing distribution per
    % the predictive distribution (derived analytically in this linear
    % gaussian case)
    p = sum(p_dyn.*repmat(p_',1,N)./repmat(normpdf(t,mso(k),sqrt(Pso(k)+q))',1,N));
    % times the filtering distribution
    p = p.*normpdf(t,mso(k),sqrt(Pso(k)));
    p = p/sum(p);
    p_ = p;
    distr(:,k) = p';
    mss(k) = t*p';
    
end
figure(1);
plot(x,mss,'-m','LineWidth',2);
pcolor(1:n,t,distr);
shading interp;
ylim([a b]);


end

if 1
% RTS smoother -- using the stationary Kalman filter

mss = zeros(1,n);
Pss = zeros(1,n);
G = Gs(20);
m = ms_st(end);
P = Ps_st(end);
mss(end) = m;
Pss(end) = P;
for k=n-1:-1:1
    m = ms_st(k)+G*(m-ms_st(k));
    P = Ps_st(k)+G^2*(P-Ps_st(k)-q);
    mss(k) = m;
    Pss(k) = P;
end
rmse(s,mss)
figure(1);
rtssolution = plot(x,mss,'--m','LineWidth',2);
hold on;
figure(2);
plot(x,Pss,'-c','LineWidth',2);

end

%% parameter estimation

% try to estimate the measurement noise
r0 = r;
N = 100;
rs = linspace(0.5,4,N);
L = zeros(1,100);

for j=1:N
    r = rs(j);
    m = 0;
    P = 1;
    m2 = 0;
    P2 = 1;
    ms = zeros(1,n);
    ms2 = zeros(1,n);
    Ps2 = zeros(1,n);
    Ps = zeros(1,n);
    Ks = Ps;
    for k=1:n
       
        % prediction
        m_= m;
        P_ = P+q;

        ms(k) = m_;
        Ps(k) = P_+r;
    end

    L(j) = exp(sum(log(normpdf(r,ms,Ps))));

end

plot(rs,L);
















