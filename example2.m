function example2
clc
n=3600; a=9; e=0.05;
ea=transpose([2+e 2 2-e 1.6+e 1.6 1.6-e 1.4+e 1.4 1.4-e 1-((a+1:n)-a)/n]);
A=diag(ea);
p=a; tau=3:8; i=8; kk=15; ll=20;
I=eye(n); X=I(:,1:p); spr=ea(1)-ea(n);
E1=zeros(ll,kk); B1a=E1; B1b=E1; C1=E1;
% data matrices for lhs and rhs in (4.6) and (4.4)
E2=zeros(ll,kk); B2a=E2; B2b=E2; C2=E2;
% data matrices for lhs and rhs in (4.7) and (4.5)
for l=1:ll
    Y=[orth(randn(p,p)); randn(n-p,p)];   % initial subspace
    Y=Y*((X'*Y)\eye(p));
    Xtau=X(:,tau); Ytau=Y(:,tau);
    Xi=X(:,1:i); Yi=Y(:,1:i);
    trighttau=sort(tan(subspacea(Xtau,Ytau)),'descend');
    trighti=sort(tan(subspacea(Xi,Yi)),'descend');
    K=orth(Y);   % initialize block Krylov subspace
    tleft=sort(tan(subspacea(Xtau,K)),'descend');
    E1(l,1)=sum(tleft);
    C1(l,1)=E1(l,1);
    psi=sort(real(eig(K'*A*K)),'descend');
    E2(l,1)=sum((ea(1:i)-psi(1:i))/spr);   % document lhs
    C2(l,1)=E2(l,1);
    sigma=zeros(p,1);
    for j=1:p   % compute sigma components with subroutine tp
        sigma(j)=1/tp(j,p,1,n,ea);
    end
    svtau=sort(sigma(tau),'descend');
    B1a(l,1)=sum(svtau(1)*trighttau);
    B1b(l,1)=sum(svtau.*trighttau);
    svi=sort(sigma(1:i),'descend');
    B2a(l,1)=sum((svi(1)*trighti).^2);
    B2b(l,1)=sum((svi.*trighti).^2);   % document rhs
    C=(2*A-(ea(p+1)+ea(n))*I)/(ea(p+1)-ea(n)); % central matrix for num_Ch
    Ka=K; Kb=C*K;   % initial subspaces for num_Ch
    for k=2:kk
        Y=orth(A*Y); K=orth([K,Y]);   % update block Krylov subspace
        if k==2
            Kc=Kb;
        else
            Kc=2*(C*Kb)-Ka;
            Ka=Kb; Kb=Kc;   % 3-term recurrence for num_Ch
        end
        tleft=sort(tan(subspacea(Xtau,K)),'descend');
        E1(l,k)=sum(tleft);
        tleft=sort(tan(subspacea(Xtau,Kc)),'descend');
        C1(l,k)=sum(tleft);
        psi=sort(real(eig(K'*A*K)),'descend');
        E2(l,k)=sum((ea(1:i)-psi(1:i))/spr);
        Kcc=orth(Kc);
        psi=sort(real(eig(Kcc'*A*Kcc)),'descend');
        C2(l,k)=sum((ea(1:i)-psi(1:i))/spr);
        for j=1:p
            sigma(j)=1/tp(j,p,k,n,ea);
        end
        svtau=sort(sigma(tau),'descend');
        B1a(l,k)=sum(svtau(1)*trighttau);
        B1b(l,k)=sum(svtau.*trighttau);
        svi=sort(sigma(1:i),'descend');
        B2a(l,k)=sum((svi(1)*trighti).^2);
        B2b(l,k)=sum((svi.*trighti).^2);   % document lhs and rhs
    end
end
e1=mean(E1); b1a=mean(B1a); b1b=mean(B1b); c1=mean(C1);
e2=mean(E2); b2a=mean(B2a); b2b=mean(B2b); c2=mean(C2);  % compute mean values
s=get(0,'ScreenSize');
figure('Position',[s(3)*0.2 s(4)*0.3 s(3)*0.6 s(4)*0.4])   % illustration
subplot(1,2,1)
semilogy(1:kk,b1a,'b - .','LineWidth',2,'MarkerSize',15)
hold on
semilogy(1:kk,b1b,'r - .','LineWidth',2,'MarkerSize',15)
semilogy(1:kk,c1,'g - .','LineWidth',2,'MarkerSize',15)
semilogy(1:kk,e1,'m - .','LineWidth',2,'MarkerSize',15)
axis([1 15 0.1^7 10^3])
xlabel('Degree of K');
legend('(4.6)','(4.4)','num Ch','num L','location','southwest')
subplot(1,2,2)
semilogy(1:kk,min(i,b2a),'b - .','LineWidth',2,'MarkerSize',15)
hold on
semilogy(1:kk,min(i,b2b),'r - .','LineWidth',2,'MarkerSize',15)
semilogy(1:kk,c2,'g - .','LineWidth',2,'MarkerSize',15)
semilogy(1:kk,e2,'m - .','LineWidth',2,'MarkerSize',15)
axis([1 15 0.1^13 10^2])
xlabel('Degree of K');
legend('(4.7)','(4.5)','num Ch','num L','location','southwest')

function s=tp(j,p,l,n,ea)   % compute Chebyshev term
d=(ea(j)-ea(p+1))/(ea(j)-ea(n));
x=(1+d)/(1-d);
te=zeros(1,l); te(1)=1; te(2)=x;
for k=3:l
    te(k)=2*x*te(k-1)-te(k-2);
end
s=te(l);