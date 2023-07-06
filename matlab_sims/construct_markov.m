function T = construct_markov(N, wbh, wwt)

Tn = [];
for i = 1:(N+1)
    x0 = i-1;
    wbar = (x0*wbh+(N-x0)*wwt)/N;
    p=(x0/N)*wbh/wbar;
    x = 0:N;
    y = binopdf(x,N,p);
    if isnan(mean(y))
            y = binopdf(x,N,1.0);
    end
    y = y';
    Tn = [Tn y];
end

T = Tn;

end