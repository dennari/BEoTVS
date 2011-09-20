%% true signal

N=200;
q=0.01^2;

s=zeros(1,N);
f=@(x)x-0.01*sin(x);
F=@(x)-0.01*cos(x);
h=@(x)0.5*sin(2*x); 
H=@(x)cos(2*x);
Y=s;
x0 = 0.4*pi;
R = [];
%x0 = 0;
x = x0;
for k=1:N
    x = f(x)+sqrt(q)*randn;
    s(k) = x;
end

% measurements
r = 0.02;
mn = 0;
n = 150;
x = floor(1:N/n:N);
y = h(s(1,x))+normrnd(0,sqrt(r),1,n);

% EKF

m = x0;
P = q;
ms = [m zeros(1,n-1)];
Ps = [P zeros(1,n-1)];
for k=2:n
    % prediction
    m_= f(m);
    P_ = F(m)^2*P+q;
    % update
    v = y(k)-h(m_);
    S = H(m_)^2*P_+r;
    K = P_*H(m_)/S;
    m = m_+K*v;
    P = P_-K^2*S;
    ms(k) = m;
    Ps(k) = P;
end
ms_ekf = ms;
Ps_ekf = Ps;
fprintf('EKF %3.4f\n',rmse(s(x),ms));
R(1) = rmse(s(x),ms);

% SLF

Ef=@(m,p)m-0.01*sin(m)*exp(-1*p/2);
Efdx=@(m,p)p-0.01*cos(m)*p*exp(-p/2);
Eh=@(m,p)0.5*sin(2*m)*exp(-2*p);
Ehdx=@(m,p)cos(2*m)*p*exp(-2*p);

m = x0;
P = q;
ms = [m zeros(1,n-1)];
Ps = [P zeros(1,n-1)];
for k=2:n
    % prediction
    m_= Ef(m,P);
    P_ = Efdx(m,P)^2/P+q;
    % update
    v = y(k)-Eh(m_,P_);
    S = Ehdx(m_,P_)^2/P_+r;
    K = Ehdx(m_,P_)/S;
    m = m_+K*v;
    P = P_-K^2*S;
    ms(k) = m;
    Ps(k) = P;
end
ms_slkf = ms;
Ps_slkf = Ps;
R(2) = rmse(s(x),ms); fprintf('SLF %3.4f\n',rmse(s(x),ms));

if 1
    plot(x,y,'.k',1:N,s,x,ms_ekf,x,ms_slkf,'MarkerSize',8);
    legend('Meas.','True','EKF','SLF');
    xlabel('t');
    ylabel('m_k');
    
    
    exportplot('ex_4_1.pdf',figW,figH,gcf,1.5);
    
    rLabels = {'RMSE'};
    cLabels = {'EKF' 'SLF'};
    matrix2latex(R,'ex_4_1_rmse.tex',...
           'alignment','d{?}{2}','format','$%.5f$','columnLabels',cLabels,...
           'rowLabels',rLabels,'rowLabelAlignment','r');
    
end
%% UKF

m = x0;
P = q;
ms = zeros(1,n);
Ps = zeros(1,n);
sig = zeros(0,3);
alpha = 12; kappa = 1;
lambda = alpha^2*(1+kappa) - 1;
Wi = 1/(2*(1+lambda));
Wm = [lambda/(1+lambda) Wi Wi];
Wc = [lambda/(1+lambda)+(1-alpha^2) Wi Wi];
for k=1:n
    % prediction
    % form the three sigma points in the original space
    sig = m*ones(3,1)+sqrt(lambda+1)*sqrt(P)*[0 1 -1]';
    sig = f(sig);
    m_ = Wm*sig;
    P_ = Wc*((sig-m_).^2)+q;
    
    sig_ = m_*ones(3,1)+sqrt(lambda+1)*sqrt(P_)*[0 1 -1]';
    sigY = h(sig_);
    
    u = Wm*sigY;
    S = Wc*((sigY-u).^2)+r;
    C = Wc*((sig_-m_).*(sigY-u));
    K = C/S;
    m = m_+K*(y(k)-u);
    P = P_-K^2*S;
    ms(k) = m;
    Ps(k) = P;
end
R(3) = rmse(s(x),ms);
if 1
    plot(1:N,s,x,ms,x,ms_slkf,'--');
    legend('True','UKF','SLF');
    xlabel('t');
    ylabel('m_k');
    exportplot('ex_5_1.pdf',figW,figH,gcf,1.5);
    rLabels = {'RMSE'};
    cLabels = {'EKF' 'SLF' 'UKF'};
    matrix2latex(R,'ex_5_1_rmse.tex',...
           'alignment','d{?}{2}','format','$%.5f$','columnLabels',cLabels,...
           'rowLabels',rLabels,'rowLabelAlignment','r');
end

%% GHKF

m = x0;
P = q;
ms = zeros(1,n);
Ps = zeros(1,n);

p = 2; % order of the polynomial
e = roots(hermite(p)); % unit sigma points
W = factorial(p)./(p*hermite(p-1,e)').^2;

for k=1:n
    % prediction
    sig = m*ones(p,1)+sqrt(P)*e;
    sig = f(sig);
    m_ = W*sig;
    P_ = W*((sig-m_).^2)+q;
    
    sig_ = m_*ones(p,1)+sqrt(P_)*e;
    sigY = h(sig_);
    
    u = W*sigY;
    S = W*((sigY-u).^2)+r;
    C = W*((sig_-m_).*(sigY-u));
    K = C/S;
    m = m_+K*(y(k)-u);
    P = P_-K^2*S;
    ms(k) = m;
    Ps(k) = P;
end
R(4) = rmse(s(x),ms);
ms_ghkf = ms;

% CKF

m = x0;
P = q;
ms = zeros(1,n);
Ps = zeros(1,n);

p = 2; % order of the polynomial
e = [1;-1]; % unit sigma points
W = ones(1,2);
W = W/sum(W);

for k=1:n
    % prediction
    sig = m*ones(p,1)+sqrt(P)*e;
    sig = f(sig);
    m_ = W*sig;
    P_ = W*((sig-m_).^2)+q;
    
    sig_ = m_*ones(p,1)+sqrt(P_)*e;
    sigY = h(sig_);
    
    u = W*sigY;
    S = W*((sigY-u).^2)+r;
    C = W*((sig_-m_).*(sigY-u));
    K = C/S;
    m = m_+K*(y(k)-u);
    P = P_-K^2*S;
    ms(k) = m;
    Ps(k) = P;
end
R(5) = rmse(s(x),ms);
ms_ckf = ms;
Ps_ckf = Ps;

if 1
    plot(1:N,s,x,ms_ckf,x,ms_ghkf,'--');
    legend('True','CKF','GHKF');
    xlabel('t');
    ylabel('m_k');
    exportplot('ex_5_2.pdf',figW,figH,gcf,1.5);
    rLabels = {'RMSE'};
    cLabels = {'EKF' 'SLF' 'UKF' 'GHKF' 'CKF'};
    matrix2latex(R,'ex_5_2_rmse.tex',...
           'alignment','d{?}{2}','format','$%.5f$','columnLabels',cLabels,...
           'rowLabels',rLabels,'rowLabelAlignment','r');
end

%% bootstrap

m = x0;
P = q;
ms = zeros(1,n);
Ps = zeros(1,n);
N = 700; % number of particles
p = m*ones(1,N); % particles
W = ones(1,N); % weights
W = W/sum(W);

for k=1:n
    % sample from the dynamic distribution, calculate weights
    %disp(W');
    for l=1:N
        p(l) = normrnd(f(p(l)),sqrt(q));
        W(l) = normpdf(y(k),h(p(l)),sqrt(r));
    end
    % normalize weights
    W = W/sum(W);
    cdf = cumsum(W);
    for l=1:N
        ran = rand;
        p(l) = p(find(cdf > ran,1));
    end
    ms(k) = W*p';
end

plot(x,ms,'-g');
fprintf('BOOTSTRAP %3.4f\n',rmse(s(x),ms));

%% SIR

m = x0;
P = q;
ms = zeros(1,n);
Ps = zeros(1,n);
N = 700; % number of particles
p = m*ones(1,N); % particles
W = ones(1,N); % weights
W = W/sum(W);

sig = zeros(0,3);
alpha = 1; kappa = 0;
lambda = alpha^2*(1+kappa) - 1;
Wi = 1/(2*(1+lambda));
Wm = [lambda/(1+lambda) Wi Wi];
Wc = [lambda/(1+lambda)+(1-alpha^2) Wi Wi];

for k=1:n
    % use UKF to get proposal mean and covariance
    sig = m*ones(3,1)+sqrt(lambda+1)*sqrt(P)*[0 1 -1]';
    sig = f(sig);
    m_ = Wm*sig;
    P_ = Wc*((sig-m_).^2)+q;
    
    sig_ = m_*ones(3,1)+sqrt(lambda+1)*sqrt(P_)*[0 1 -1]';
    sigY = h(sig_);
    
    u = Wm*sigY;
    S = Wc*((sigY-u).^2)+r;
    C = Wc*((sig_-m_).*(sigY-u));
    K = C/S;
    m = m_+K*(y(k)-u);
    P = P_-K^2*S;
   
    % sample from the proposal distribution
    for l=1:N
        pp = normrnd(m,P);
        W(l) = W(l)*normpdf(y(k),h(pp),sqrt(r))*normpdf(pp,f(p(l)),sqrt(q))/normpdf(pp,m,P);
        p(l) = pp;
    end
    % normalize weights
    W = W/sum(W);
    % resample if needed
    Neff = 1/(W*W');
    if Neff < N / 10
        %fprintf('resampling, neff = %3.1f\n',Neff);
        po = p;
        cdf = cumsum(W);
        for l=1:N
            ran = rand;
            p(l) = po(find(cdf > ran,1));
        end
        W = ones(1,N); % weights
        W = W/sum(W);
    else
        %fprintf('not resampling, neff = %3.1f\n',Neff);
    end
    ms(k) = W*p';
    

end

plot(x,ms,'-r');
fprintf('SIR %3.4f\n\n',rmse(s(x),ms));

legend('Meas.','True','Bootstrap','SIR');
xlabel('t');
ylabel('m_k');
fprintf('SLF %3.4f\n',rmse(s(x),ms));
exportplot('ex_6_2.pdf',figW,figH);


%% ERTS

mss = zeros(1,n);
Pss = mss;
m_s = ms_ekf(:,end);
P_s = Ps_ekf(:,end);
mss(:,end) = m_s;
Pss(:,end) = P_s;

for k=n-1:-1:1
    m = ms_ekf(k);
    P = Ps_ekf(k);
    % prediction
    m_ = f(m);
    P_ = F(m)^2*P+q;
    
    % update
    G = P*F(m)/P_;
    m_s = m + G*(m_s-m_);
    P_s = P + G*(P_s-P_)*G';
    mss(k) = m_s;
    Pss(k) = P_s;
    
end

plot(x,mss,'--r','LineWidth',2);
fprintf('ERTS %3.4f\n\n',rmse(s(x),mss));

%% SL RTS

mss = zeros(1,n);
Pss = mss;
m_s = ms_slkf(:,end);
P_s = Ps_slkf(:,end);
mss(:,end) = m_s;
Pss(:,end) = P_s;

for k=n-1:-1:1
    m = ms_slkf(k);
    P = Ps_slkf(k);
    % prediction
    m_= Ef(m,P);
    P_ = Efdx(m,P)^2/P+q;
    
    % update
    G = Efdx(m,P)/P_;
    m_s = m + G*(m_s-m_);
    P_s = P + G*(P_s-P_)*G';
    mss(k) = m_s;
    Pss(k) = P_s;
       
end

plot(x,mss,'--g','LineWidth',2);
fprintf('SLRTS %3.4f\n\n',rmse(s(x),mss));

legend('Meas.','True','EKF','SLF','ERTS','SL RTS');
xlabel('t');
ylabel('m_k');
fprintf('SLF %3.4f\n',rmse(s(x),ms));
exportplot('ex_7_3.pdf',figW,figH);


