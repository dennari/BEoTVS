%%
syms N r w q b111 b112 b121 b122 b211 b212 b221 b222 b311 b312 b321 b322 real;

B1 = [b111 b112; b121 b122];
B2 = [b211 b212; b221 b222];
B3 = [b311 b312; b321 b322];


Q=[(q*w-q*cos(w)*sin(w))/(2*w^3) (q*sin(w)^2)/(2*w^2);...
   (q*sin(w)^2)/(2*w^2) (q*w+q*cos(w)*sin(w))/(2*w)];
A=[cos(w) sin(w)/w; -w*sin(w) cos(w)];

w = 0.5
q = 0.01


L = -N/2*log(det(Q))-0.5*trace(Q^-1*(B1-B2*A'-A*B2+A*B3*A'));
Lw = subs(subs(subs(subs(subs(L,B1,ones(2,2)),B2,ones(2,2)),B3,ones(2,2)),q,0.01),N,100);
Lq = subs(subs(subs(subs(subs(L,B1,ones(2,2)),B2,ones(2,2)),B3,ones(2,2)),w,10),N,100);


%%
x = linspace(-10,10,5000);
y = subs(diff(Lw),x);
plot(x,y)
y = subs(diff(Lq),x);
figure;
plot(x,y)
