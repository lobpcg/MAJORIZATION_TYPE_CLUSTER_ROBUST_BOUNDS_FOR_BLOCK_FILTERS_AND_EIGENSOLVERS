function thm3_5
clc
for l=1:1000
    n=100; majleTol=n*n*eps;
    X=orth(randn(n)+1i*randn(n));   % eigenbasis
    % ea=sort(randn(n,1),'descend');
    ef=randn(n,1)+1i*randn(n,1);   % f(lambda)
    p=10;
    r=max(abs(ef(p+1:n)))/min(abs(ef(1:p)));
    ef(1:p)=(r*1.01)*ef(1:p);   % assumption (2.1)
    % A=X*diag(ea)*X';
    F=X*diag(ef)*X';   % f(A)
    tau=unique(ceil(rand(1,9)*p));
    X=X(:,1:p);
    Yt=orth(randn(n,p)+1i*randn(n,p));
    Y=Yt*((X'*Yt)\eye(p));   % Lemma 3.1
    Yp=F*Y;
    Xtau=X(:,tau);
    Ytau=Y(:,tau);
    tleft=sort(tan(subspacea(Xtau,Yp)),'descend');
    tright=sort(tan(subspacea(Xtau,Ytau)),'descend');
    phitau=sort(1./abs(ef(tau)),'descend');
    phih=sort(abs(ef(p+1:n)),'descend');
    phiht=phih(1:length(tau));
    ifm=majle(tleft,phitau.*phiht.*tright,majleTol);
    if ifm==0
        disp([tleft, phitau.*phiht.*tright])
    end
end