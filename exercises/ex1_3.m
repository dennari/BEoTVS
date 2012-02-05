%% true signal

A = [1 1;0 1];
Q = diag([1/100 1]);
H = [1 0];
R = 100;

n = 100;
x0 = mvnrnd([0;0],eye(2))';

X = [x0 zeros(2,n-1)];
for k=2:n
    X(:,k) = A*X(:,k-1)+mvnrnd([0;0],Q)';
end
Y = X(1,:)+sqrt(R)*randn(1,n);

%% Filter

MM = zeros(2,n);
PP = zeros(2,2,n);
M = zeros(2,1);
P = eye(2);
for k=1:size(Y,2)
    %
    % Track with KF
    %
    [M,P] = kf_predict(M,P,A,Q);
    [M,P] = kf_update(M,P,Y(k),H,R);

    MM(:,k) = M;
    PP(:,:,k) = P;
end

if 1 %plot
    plot(1:n,X(1,:),1:n,Y,'k.',1:n,MM(1,:),'--','MarkerSize',8);
    set(gca,'YTick',[]);
    legend('Signal 1st','Measurements','Filter 1st mean');
    xlabel('Time');
    exportplot('ex_1_3_signal.pdf',figW,figH,gcf,1.5);
    figure;
    plot(1:n,X(2,:),1:n,MM(2,:),'--');
    set(gca,'YTick',[]);
    legend('Signal 2nd','Filter 2nd mean');
    xlabel('Time');
    exportplot('ex_1_3_derivative.pdf',figW,figH,gcf,1.5);
    
    rLabels = {'RMSE'};
    cLabels = {'Filter' 'Measurements'};
    matrix2latex([rmse(X(1,:),MM(1,:)) rmse(X(1,:),Y)],'ex_1_3_rmse.tex',...
           'alignment','d{?}{2}','format','$%.5f$','columnLabels',cLabels,...
           'rowLabels',rLabels,'rowLabelAlignment','r');
end
