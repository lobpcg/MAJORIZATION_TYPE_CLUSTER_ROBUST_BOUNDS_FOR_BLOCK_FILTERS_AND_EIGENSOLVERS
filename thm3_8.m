function thm3_8
clc
for l=1:1000
    n=100; majleTol=n*n*eps;
    X=orth(randn(n)+1i*randn(n));   % eigenbasis
    ea=sort(randn(n,1),'descend');   % lambda
    ef=randn(n,1)+1i*randn(n,1);   % f(lambda)
    p=10;
    r=max(abs(ef(p+1:n)))/min(abs(ef(1:p)));
    ef(1:p)=(r*1.01)*ef(1:p);   % assumption (2.1)
    A=X*diag(ea)*X';
    F=X*diag(ef)*X';   % f(A)
    i=6; ii=1:i;
    X=X(:,1:p);
    Yt=orth(randn(n,p)+1i*randn(n,p));
    Y=Yt*((X'*Yt)\eye(p));   % Lemma 3.1
    Yp=F*Y;
    V=orth(Yp); eta=sort(real(eig(V'*A*V)),'descend');
    Xi=X(:,ii);
    Yi=Y(:,ii);
    delta=sort((ea(ii)-eta(ii))./(eta(ii)-ea(n)),'descend');
    tright1=sort(tan(subspacea(Xi,Yi)),'descend');
    phii=sort(1./abs(ef(ii)),'descend');
    phih=sort(abs(ef(p+1:n)),'descend');
    phihi=phih(ii);
    ifm=majle(delta,(phii.*phihi.*tright1).^2,majleTol);
    if ifm==0
        disp([delta,(phii.*phihi.*tright1).^2])
    end
    tright=sort(tan(subspacea(X,Y)),'descend');
    tright2=tright(ii);
    ifm=majle(delta,(phii.*phihi.*tright2).^2,majleTol);
    if ifm==0
        disp([delta,(phii.*phihi.*tright2).^2])
    end
end