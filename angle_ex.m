function []=angle_ex(figW,figH)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Bayesian Estimation of Time-Varying Processes (5 p) L V
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Generate data
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    randn('state',123);

    %
    % Create a bit curved trajectory and angle
    % measurements from two sensors
    %

    S = [-1.5 1;0.5 1];


    sd = 0.05;      % Standard deviation of measurements
    dt = 0.01;      % Sampling period
    x0 = [0;0;1;0]; % Initial state
    len = 3;
    a = zeros(1,500);
    a(1,50:100)  = pi/2/51/dt + 0.01*randn(1,51);
    a(1,200:250) = pi/2/51/dt + 0.01*randn(1,51);
    a(1,350:400) = pi/2/51/dt + 0.01*randn(1,51);
    x = x0;
    t = 0;
    X = [];
    Theta = [];
    T = [];
    %% generate measurements
    for i=1:500
        F = [0 0  1    0;...
             0 0  0    1;...
             0 0  0   a(i);...
             0 0 -a(i) 0];
        x = expm(F*dt)*x;
        y = h(x)+normrnd(0,sd,size(S,2),1); 
        t  = t + dt;
        X = [X x];
        T = [T t];
        Theta = [Theta y];
    end
    steps = size(Theta,2);

    %% common variables
    m = [];
    P = [];
    k = [];
    ms = zeros(4,steps);
    Ps = zeros(4,4,steps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Parameters of the dynamic model
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qc = 0.1;

    %% This is the transition matrix
    A  = [1 0 dt 0;
        0 1 0 dt;
        0 0 1 0;
        0 0 0 1];

    %% This is the process noise covariance
    Q = [qc*dt^3/3 0 qc*dt^2/2 0;
       0 qc*dt^3/3 0 qc*dt^2/2;
       qc*dt^2/2 0 qc*dt 0;
       0 qc*dt^2/2 0 qc*dt];



    R  = sd^2*eye(size(S,2));   % The joint covariance
    rLabels = {'RMSE'};
    if 1
        f1 = figure;ax(1) = gca;
        M1 = baseline();
        axis equal;
        legend('True','Baseline');


        f2 = figure;ax(2) = gca;
        linkaxes(ax);
        M2 = ekf();
        axis equal;
        legend('True','EKF');    
        xlim([-2 2]);ylim([-2.2 2]);

        exportplot('ex_4_3_baseline.pdf',figW,figH,f1,1.5);
        exportplot('ex_4_3_ekf.pdf',figW,figH,f2,1.5);
           
        cLabels = {'Baseline' 'EKF'};
        matrix2latex([rmse(X,M1) rmse(X,M2)],'ex_4_3_rmse.tex',...
               'alignment','d{?}{2}','format','$%.5f$','columnLabels',cLabels,...
               'rowLabels',rLabels,'rowLabelAlignment','r');
    end
    if 0
        f3 = figure;ax(1) = gca;
        M3 = ckf();
        axis equal;
        legend('True','CKF/UKF');

        f4 = figure;ax(2) = gca;
        linkaxes(ax);
        M4 = ukf();
        axis equal;
        legend('True','UKF');    
        xlim([-2 2]);ylim([-2.2 2]);

        exportplot('ex_5_3_ckfukf.pdf',figW,figH,f3);
        exportplot('ex_5_3_ukf.pdf',figW,figH,f4);
        
        cLabels = {'CKF' 'UKF'};
        matrix2latex([rmse(X,M3) rmse(X,M4)],'ex_5_3_rmse.tex',...
               'alignment','d{?}{2}','format','$%.5f$','columnLabels',cLabels,...
               'rowLabels',rLabels,'rowLabelAlignment','r');
    end
    if 0
        f5 = figure;ax(1) = gca;
        bootstrap();
        axis equal;
        legend('True','Bootstrap');

        f6 = figure;ax(2) = gca;
        linkaxes(ax);
        sir();
        axis equal;
        legend('True','SIR');    
        xlim([-2 2]);ylim([-2.2 2]);

        exportplot('ex_6_3_bootstrap.pdf',figW,figH,f5);
        exportplot('ex_6_3_sir.pdf',figW,figH,f6);
    end
   
    function [o,h]=baseline()
        m = x0;     % Initialize to true value
        ms = zeros(4,steps);

        for k=1:steps
            %% Compute crossing of the measurements
            dx1 = cos(Theta(1,k));
            dy1 = sin(Theta(1,k));
            dx2 = cos(Theta(2,k));
            dy2 = sin(Theta(2,k));
            d = [dx1 dx2; dy1 dy2]\[S(1,2)-S(1,1);S(2,2)-S(2,1)];
            cross_xy = S(:,1) + [dx1;dy1]*d(1);
            %% Compute estimate
            m(3:4) = [0;0];
            m(1:2) = cross_xy;
            ms(:,k) = m;
        end

        fprintf('BL %3.4f\n',rmse(X,ms));
        h= showtrace('b');
        %disp(h);
        o = ms;
    end

    %%%% EKF %%%%%%%%%%
    
    function [o,hh]=ekf()
        m = x0;            % Initialize to true value
        P = eye(4);        % Some uncertainty
        ms = zeros(4,steps);
        for k=1:steps
            %% Compute estimate here
            m_ = A*m;
            P_ = A*P*A'+ Q;
            y = Theta(:,k); 
            v =  y - h(m_);
            S_ = H(m_)*P_*H(m_)'+R;
            K = P_*H(m_)'/S_;
            m = m_+K*v;
            P = P_-K*S_*K';
            ms(:,k) = m;
            Ps(:,:,k) = P;
        end

        %% Compute error
        fprintf('EKF %3.4f\n',rmse(X,ms));
        o = ms;
        hh = showtrace('b');
        axis([-2 2 -2.5 1.5]);
    end

    function o=h(x)
        o = atan2(x(2)-S(2,:)',x(1)-S(1,:)');
    end

    function o=H(x)
        H1 = @(xx,s)-(xx(2)-s(2))/((xx-s)'*(xx-s));
        H2 = @(xx,s)(xx(1)-s(1))/((xx-s)'*(xx-s));
        o = zeros(size(S,2),4);
        for ii = 1:size(S,2)
            o(ii,1:2) = [H1(x(1:2),S(:,ii)) H2(x(1:2),S(:,ii))];
        end    
    end



    function ms = ckf()

        m = x0;            % Initialize to true value
        P = eye(4);        % Some uncertainty
        ms = zeros(4,steps);

        for k=1:steps
            [m,P] = ckf_(Theta(:,k));
            ms(:,k) = m;
            Ps(:,:,k) = P;
        end

          %% Compute error
        fprintf('CKF %3.4f\n',rmse(X,ms));
        showtrace('b');
    end

    function [mo,Po] = ckf_(yy)
        n = 4;
        nn = size(S,2);
        e = [sqrt(n)*eye(n) -sqrt(n)*eye(n)]; % unit sigma points
        W = ones(1,2*n);
        W = W/sum(W);
        
        sig = repmat(m,1,2*n)+chol(P,'lower')*e;
        sig = A*sig;
        m_ = sig*W';
        m__ = repmat(m_,1,2*n);
        P_ = (sig-m__)*(sig-m__)'/(2*n)+Q;
        sig_ = m__+chol(P_,'lower')*e;
        sigY = zeros(nn,2*n);
        for j=1:2*n
            sigY(:,j)=h(sig(:,j));
        end
        u = sigY*W';
        u__ = repmat(u,1,2*n);
        S_ = (sigY-u__)*(sigY-u__)'/(2*n)+R;
        C = (sig_-m__)*(sigY-u__)'/(2*n);
        K = C/S_;
        mo = m_+K*(yy-u);
        Po = P_-K*S_*K';
    end




    function ms=ukf()
        m = x0;            % Initialize to true value
        P = eye(4);        % Some uncertainty
        ms = zeros(4,steps);
        n = 4;
        nn = size(S,2);
        alpha = 1; kappa = 0;
        lambda = alpha^2*(1+kappa) - 1;
        e = [zeros(n,1) sqrt(n+lambda)*eye(n) -sqrt(n+lambda)*eye(n)]; % unit sigma points
        Wm = [lambda/(n+lambda) 1/(2*(n+lambda))*ones(1,2*n)];
        Wc = [lambda/(n+lambda)+(1-alpha^2) 1/(2*(n+lambda))*ones(1,2*n)];
        for k=1:steps
            % prediction
            sig = repmat(m,1,2*n+1)+chol(P,'lower')*e;
            sig = A*sig;
            m_ = sig*Wm';
            m__ = repmat(m_,1,2*n+1);
            P_ = (sig-m__)*diag(Wc)*(sig-m__)'+Q;
            sig_ = m__+chol(P_,'lower')*e;
            sigY = zeros(nn,2*n);
            for j=1:(2*n+1)
                sigY(:,j)=h(sig(:,j));
            end
            u = sigY*Wm';
            u__ = repmat(u,1,2*n+1);
            S_ = (sigY-u__)*diag(Wc)*(sigY-u__)'+R;
            C = (sig_-m__)*diag(Wc)*(sigY-u__)';
            K = C/S_;
            m = m_+K*(Theta(:,k)-u);
            P = P_-K*S_*K';
            ms(:,k) = m;
            Ps(:,:,k) = P;
        end

        fprintf('UKF %3.4f\n',rmse(X,ms));
        showtrace();
    end

    function ms=bootstrap()
        P = eye(4);        % Some uncertainty
        ms = zeros(4,steps);
        ii = zeros(1,steps);
        N = 150;
        p = repmat(x0,1,N); % particles
        W = ones(1,N); % weights
        W = W/sum(W);

        for k=1:steps
            % prediction
            for l=1:N
                p(:,l) = mvnrnd(A*p(:,l),Q);
                W(l) = mvnpdf(Theta(:,k),h(p(:,l)),R);
            end
            W = W/sum(W);
            cdf = cumsum(W);
            po = p;
            for l=1:N
                ran = rand;
                ii(k) = find(cdf >= ran,1);
                p(:,l) = po(:,find(cdf >= ran,1));
            end
            ms(:,k) = p*W';
            %anim();
        end

        %% Compute error
        fprintf('BOOTSTRAP %3.4f\n',rmse(X,ms));
        showtrace();
    end

    function ms=sir()
        m = x0;
        P = eye(4);        % Some uncertainty
        ms = zeros(4,steps);
        ii = zeros(1,steps);
        N = 150;
        p = repmat(x0,1,N); % particles
        W = ones(1,N); % weights
        W = W/sum(W);
        c = 0;
        for k=1:steps
            
            [m,P] = ckf_(Theta(:,k));
            % prediction
            for l=1:N
                pn = mvnrnd(m,P)';
                po = p(:,l);
                W(l) = W(l)*mvnpdf(Theta(:,k),h(pn),R)*mvnpdf(pn,A*po,Q)/mvnpdf(pn,m,P);
                if(~W(l))
                    W(l) = 1/N;
                    c=c+1;
                end
                p(:,l) = pn;
            end
            W_ = W/sum(W);
            if sum(isnan(W_)) > 0
                Theta(:,k)
                h(pn)
                mvnpdf(Theta(:,k),h(pn),R)
                pn
                A*po
                mvnpdf(pn,A*po,Q)
                mvnpdf(pn,m,P)
                W
                k
                break;
            end
            W = W_;
            % resample if needed
            Neff = 1/(W*W');
            if Neff < N / 10
                %fprintf('resampling, neff = %3.1f\n',Neff);
                po = p;
                cdf = cumsum(W);
                for l=1:N
                    ran = rand;
                    p(:,l) = po(:,find(cdf > ran,1));
                end
                W = ones(1,N); % weights
                W = W/sum(W);
            else
                %fprintf('not resampling, neff = %3.1f\n',Neff);
            end
            ms(:,k) = p*W';
            %anim();
        end
        c
        %% Compute error
        fprintf('SIR %3.4f\n',rmse(X,ms));
        showtrace();
    end


    function []=diff()
      Hest1 = zeros(1,size(X,2));
      Hest2 = Hest1;
      for k=2:30:size(X,2)
          y = Theta(:,k);
          y_ = Theta(:,k-1);
          x = X(:,k);
          x_ = X(:,k-1);
          %(x(1)-x_(1))
          %H(x)'
          %[(y(1)-y_(1))/(x(1)-x_(1)) (y(1)-y_(1))/(x(2)-x_(2)) 0 0;
          % (y(2)-y_(2))/(x(1)-x_(1)) (y(2)-y_(2))/(x(2)-x_(2)) 0 0]'   

          Hest1(k) = [1 0]*H(x_)*[1 0 0 0]';
          Hest2(k) = (y(1)-y_(1))/(x(1)-x_(1));

      end
      subplot(3,1,1);
      plot(Hest1);
      subplot(3,1,2);
      Hest2(isinf(Hest2)) = 0;
      plot(Hest2);
      subplot(3,1,3);
      plot(X(1,:),Theta(1,:));
    end
  


    function []=anim()
        if rem(k,10) == 1
            clf;
            plot(X(1,:),X(2,:),'r-',...
               m(1),m(2),'bo',...
               ms(1,1:k),ms(2,1:k),'b--');
            hold on; 
            for ii=1:size(S,2)
              plot([S(1,ii);S(1,ii)+len*cos(Theta(ii,k))],[S(2,ii);S(2,ii)+len*sin(Theta(ii,k))],'k--');
            end
            hold off;
            axis([-2 2 -2.5 1.5]);
            drawnow;
        end
    end
    function [h]=showtrace(color)
        %clf;
        if ~nargin
            color='b';
        end
%         l95 = zeros(2,steps);
%         u95 = l95;
%         l95(1,:) = ms(1,:)-2*sqrt(squeeze(Ps(1,1,:))');
%         l95(2,:) = ms(2,:)-2*sqrt(squeeze(Ps(2,2,:))');
%         u95(1,:) = ms(1,:)+2*sqrt(squeeze(Ps(1,1,:))');
%         u95(2,:) = ms(2,:)+2*sqrt(squeeze(Ps(2,2,:))');
%         
%         
%        h = plot(X(1,:),X(2,:),'r-',...
%                 ms(1,:),ms(2,:),[color '-'],...
%                 l95(1,:),l95(2,:),'--k',...
%                 u95(1,:),u95(2,:),'--k',...
%                 S(1,:)',S(2,:)','kx');
            
       h = plot(X(1,:),X(2,:),'r-',...
                ms(1,:),ms(2,:),[color '-'],...
                S(1,:)',S(2,:)','kx');     
            
       %disp(h);
       axis([-2 2 -2.5 1.5]);
       
    end
end



  

  
  
  
  

