function o = energy(w,q,r,y)
	s = sin(w);
	c = cos(w);
	Q = [(q*w-q*c*s)/(2*w^3), q*s^2/(2*w^2);
		  q*s^2/(2*w^2),	 (q*w+q*c*s)/(2*w)];
	A = [c s/w;
		-w*s c];
	H = [0 1];
	xDim = size(A,1);
	yDim = size(y,1);
	% we suppose the data is given as column vectors, i.e y \in R^(2xN)
	N = size(y,2);
	
	%ms = zeros(xDim,N);
	%Ps = zeros(xDim,xDim,N);
	
	m = zeros(xDim,1);
	P = eye(xDim);


	%ms(:,1) = m;
	%Ps(:,:,1) = P;
	o = 0;
	for k=1:N
		m_ = A*m;
		P_ = A*P*A'+Q;
		S = H*P_*H' + r;
		K = P_*H'/S;

		% p(x_k|y_1:k)=N(m,P)
		d = y(:,k)-H*m_;
		m = m_+K*d;
		P = P_-K*S*K';

		% p(y_k|y_1:k-1)=N(H*m_,S)
		%ms(:,k) = H*m_;
		%Ps(:,:,k) = S;

		o = o + log(det(2*pi*S))+d^2/S;
	end 
	o = 0.5*o;

	


end