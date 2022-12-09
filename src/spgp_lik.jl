function spgp_lik(w,y,x,n,del=1e-6)
    [N,dim] = size(x);
    xb = reshape(w[1:end-dim-2],n,dim);
    b = exp(w[end-dim-1:end-2]);
    c = exp(w[end-1]);
    sig = exp(w[end]);

    xb = xb.*repeat(transpose(sqrt(b)),n,1);
    x = x.*repeat(transpose(sqrt(b)),N,1);

    Q = xb*transpose(xb);
    Q = repeat(diagm(Q),1,n) + repeat(diagm(transpose(Q), n, 1) - 2*Q);
    Q = c*exp(-0.5*Q) + del*I(n);

    K = -2*xb*transpose(x) + repeat(transpose(sum(x.*x,2)), n, 1) + repeat(sum(xb.*xb,2), 1, N);
    K = c*exp(-0.5*K);

    L = transpose(cholesky(Q));
    V = L\K;
    ep = 1 + (c-transpose(sum(V.^2)))/sig;
    K = K./repeat(transpose(sqrt(ep)), n, 1);
    V = V./repeat(transpose(sqrt(ep)), n, 1);
    y = y./sqrt(ep);
    Lm = transpose(cholesky(sig*I(n) + V*transpose(V)));
    invLmV = Lm\V;
    bet = invLmV*y;

    fw = sum(log(diagm(Lm))) + (N-n)/2*log(sig) + (
        transpose(y)*y - transpose(bet)*bet)/2/sig + sum(log(ep))/2 + 0.5*N*log(2*pi); 
    
    Lt = L*Lm;
    B1 = transpose(Lt)\invLmV;
    b1 = transpose(Lt)\bet;
    invLv = transpose(L)\V;
    let invL = inv(L)
    invQ = transpose(invL)*invL;
    end
    let invLt = inv(Lt)
    invA = transpose(invLt)*invLt;
    end
    mu = transpose(transpose(transpose(Lm)\bet)*V);
    sumVsq = transpose(sum(V.^2));
    bigsum = y.*(bet'*invLmV)'/sig - sum(invLmV.*invLmV)'/2 - (y.^2+mu.^2)/2/sig + 0.5;
    TT = invLV*(invLV'.*repeat(bigsum,1,n));

    for i = 1:dim
        dnnQ = dist(xb[:,i],xb[:,i]).*Q;
        dNnK = dist(-xb[:,i],-x[:,i]).*K;

        dfxb[:,i] = - b1.*(dNnK*(y-mu)/sig + dnnQ*b1) + sum((invQ - invA*sig).*dnnQ,2) + epdot*bigsum - 2/sig*sum(dnnQ.*TT,2); 
        dfb[i,1] = (((y-mu)'.*(b1'*dNnK))/sig + (epPmod.*bigsum)')*x[:,i];
        dNnK = dNnK.*B1; 
        dfxb[:,i] = dfxb[:,i] + sum(dNnK,2);
        dfb[i,1] = dfb[i,1] - sum(dNnK,1)*x[:,i];
        dfxb[:,i] = dfxb[:,i]*sqrt(b[i]);
        dfb[i,1] = dfb[i,1]/sqrt(b[i]);
        dfb[i,1] = dfb[i,1] + dfxb[:,i]'*xb[:,i]/b[i];
        dfb[i,1] = dfb[i,1]*sqrt(b[i])/2;
    end

    epc = (c./ep - sumVsq - del*sum((invLV).^2)')/sig;
    
    dfc = (n + del*trace(invQ-sig*invA) - sig*sum(sum(invA.*Q')))/2 - mu'*(y-mu)/sig + b1'*(Q-del*eye(n))*b1/2 + epc'*bigsum;
    dfsig = sum(bigsum./ep);
    
    dfw = [reshape(dfxb,n*dim,1);dfb;dfc;dfsig];

    return fw, dfw
end

