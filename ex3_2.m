% for outputting figures correctly
latexTextW = 341; % in points
figAspectRatio = 16/10; % w/h
figW = 1.3*latexTextW;
figH = figW/figAspectRatio;

%% Exercise round 3, exercise 2

%true signal
N=100;
q=3;

% signal
s=cumsum([0 normrnd(0,sqrt(q),1,N-1)]);

% measurements
A = 1;
H = 1;
r = 5;
mn = 0;
n = N;
x = floor(1:N/n:N);
y = s(1,x)+normrnd(mn,sqrt(r),1,n);

%% Kalman filter
m = 0;
P = 1;
ms = zeros(1,n);ms(1)=m;
Ps = zeros(1,n);Ps(1)=P;
Ks = zeros(1,n-1);
for k=2:n
    % update
    K = (P+q)/(P+q+r);
    m = (1-K)*m+K*y(k);
    P = P+q-K^2*(P+q+r);
    ms(k) = m;
    Ps(k) = P;
    Ks(k-1) = K;
end

mso = ms;
Pso = Ps;
rmse(x,ms)

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
dms = ms;
dPs = Ps; 

if 1 % plot
    % means, measurements
    plot(x,y,'.k',...
        1:n,s,...
        x,dms,...
        x,mso);
    legend('Measurements','True','Discretization','Kalman');
    xlabel('Time');
    ylabel('Signal mean');
    exportplot('ex_3_2_means.pdf',figW,figH,gcf,1.5);
    % variances
    figure;
    plot(x,dPs,x,Pso,'--');
    legend('Kalman','Discretization');
    xlabel('Time');
    ylabel('Variance');
    exportplot('ex_3_2_variances.pdf',figW,figH,gcf,1.5);
end


%% stationary Kalman filter (use the last value of K (the Kalman gain) from
% the normal filter)

m = 0;
P = 1;
ms = zeros(1,n);ms(1)=m;
Ps = zeros(1,n);Ps(1)=P;
for k=2:n
    m = (1-K)*m+K*y(k);
    P = P+q-K^2*(P+q+r);
    ms(k) = m;
    Ps(k) = P;
end
%kalmansolution = plot(x,ms,'-k');
ms_st = ms;
Ps_st = Ps;
rmse(x,ms)




%% RTS smoother

mss = zeros(1,n);
Pss = zeros(1,n);
Gs = mss(1:n-1);
m = mso(end); % mso = Kalman filter means
P = Pso(end); % Pso = Kalman filter variances
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
rmse(s,mss)
rts_m = mss;
rts_p = Pss;

if 1 % plot
    % compare with kalman filter
    % means, measurements
    plot(x,y,'.k',...
        1:n,s,...
        x,mso,...
        x,mss);
    legend('Measurements','True','Kalman','RTS');
    xlabel('Time');
    ylabel('Signal mean');
    exportplot('ex_7_1a_means.pdf',figW,figH,gcf,1.5);

    % variances
    figure;
    plot(x,Pso,x,Pss,'--');
    legend('Kalman','RTS');
    xlabel('Time');
    ylabel('Variance');
    exportplot('ex_7_1a_variances.pdf',figW,figH,gcf,1.5);
    
end


if 1
% RTS smoother -- discrete

mss = zeros(1,n);
Pss = zeros(1,n);
m = mso(end);
P = Pso(end);
mss(end) = m;
Pss(end) = P;
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
    m = t*p';
    mss(k) = m;
    Pss(k) = (t-m).^2*p';
end
rts_discr_m = mss;
rts_discr_P = Pss;

if 1 %plot
    %pcolor(1:n,t,distr);
    %hold on;
    %shading interp;
    %ylim([a b]);
    plot(x,rts_discr_m,x,rts_m,'--');
    legend('Discretization','RTS');
    xlabel('Time');
    ylabel('Signal mean');
    exportplot('ex_7_1b_means.pdf',figW,figH,gcf,1.5);
    figure;
    plot(x,rts_discr_P,x,rts_p,'--');
    legend('Discretization','RTS');
    xlabel('Time');
    ylabel('Variance');
    exportplot('ex_7_1b_variances.pdf',figW,figH,gcf,1.5);
end

end

%%

if 1
% RTS smoother -- using the stationary Kalman filter

mss = zeros(1,n);
Pss = zeros(1,n);
G = Gs(end);
m = ms_st(end); % ms_st = stationary Kalman filter means
P = Ps_st(end); % Ps_st = stationary Kalman filter variances
mss(end) = m;
Pss(end) = P;
for k=n-1:-1:1
    m = ms_st(k)+G*(m-ms_st(k));
    P = Ps_st(k)+G^2*(P-Ps_st(k)-q);
    mss(k) = m;
    Pss(k) = P;
end
rmse(s,mss)
if 1 %plot
    plot(x,mss,x,rts_m,'--');
    legend('Stat. RTS','RTS');
    xlabel('Time');
    ylabel('Signal mean');
    exportplot('ex_7_1c_means.pdf',figW,figH,gcf,1.5);
    figure;
    plot(x,Pss,x,rts_p,'--');
    legend('Stat. RTS','RTS');
    xlabel('Time');
    ylabel('Variance');
    exportplot('ex_7_1c_variances.pdf',figW,figH,gcf,1.5);
end

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
















