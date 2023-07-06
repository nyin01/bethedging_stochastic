function T = p1_meanmat(wa, wb, s, pSpec, N, emodel)
gmf = geomean([wa, wb]) + s;
wc = (5*pSpec - 2*(4*gmf^2 + (9*pSpec^2)/4)^(1/2))/(4*pSpec - 4);

x = 0:N;
Ta = [];
Tb = [];
yy0 = zeros(N+1, 1);
yy0(1) = 1;
Ta = [Ta yy0];
Tb = [Tb yy0];
for i = 1:N
    x0 = i;
    allc = [0:x0]; %all options for realized bh counts with each phenotype
    phenodist = binopdf(allc, x0, pSpec);
    r_pS = allc./x0; % all possibilities for realized pS
    
    tmata = [];
    tmatb = [];
    for j = 1:length(r_pS)
        pS = r_pS(j);
        
        %for the p1 model
        %env A
        wbh_a = wc*(1-pS)+wa*(pS);
        wbar_a = (x0*wbh_a+(N-x0)*wa)/N;
        p_a=(x0/N)*wbh_a/wbar_a;
        y_a = binopdf(x,N,p_a);
        if isnan(mean(y_a))
            y_a = binopdf(x,N,1.0);
        end
        y = y_a.*phenodist(j);
        tmata = [tmata; y];
        
        % env b
        wbh_b = wc*(1-pS)+wb*(pS);
        wbar_b = (x0*wbh_b+(N-x0)*wb)/N;
        p_b=(x0/N)*wbh_b/wbar_b;
        y_b = binopdf(x,N,p_b);
        if isnan(mean(y_b))
            y_b = binopdf(x,N,1.0);
        end
        yb = y_b.*phenodist(j);
        tmatb = [tmatb; yb];
        
        
    end
    
    yya = sum(tmata);
    Ta = [Ta yya'];
    
    yyb = sum(tmatb);
    Tb = [Tb yyb'];
end

if emodel == 0
    AB = Tb*Ta;
    BA = Ta*Tb;
    T = (AB+BA)./2;
elseif emodel == 1
    AB = Tb*Ta;
    BA = Ta*Tb;
    AA = Ta*Ta;
    BB = Tb*Tb;
    T = (AB+BA+AA+BB)./4;
end
    

end