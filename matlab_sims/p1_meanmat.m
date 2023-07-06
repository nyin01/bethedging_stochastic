function T = p1_meanmat(wa, wb, s, pSpec, N, emodel)
%function for both the E1P1 and E0P1 models
%final argument is whether env stoch is on 1 or off 0

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
    x0 = i; %phenotypic distribution is a function of the number of bh individuals in the population
    allc = [0:x0]; %all options for realized bh counts with each phenotype
    phenodist = binopdf(allc, x0, pSpec); %probability of experiencing every possible realized pS
    r_pS = allc./x0; % all possibilities for realized pS
    
    tmata = [];
    tmatb = [];
    for j = 1:length(r_pS) %for every possible realized pS --> what is the fitness of the BH --> what is the transition vector
        pS = r_pS(j);
        
        %env A
        wbh_a = wc*(1-pS)+wa*(pS);
        wbar_a = (x0*wbh_a+(N-x0)*wa)/N;
        p_a=(x0/N)*wbh_a/wbar_a;
        y_a = binopdf(x,N,p_a);
        if isnan(mean(y_a))
            y_a = binopdf(x,N,1.0);
        end
        y = y_a.*phenodist(j); %weight transition probability vector by probability of having this realized pS
        tmata = [tmata; y];
        
        % env b
        wbh_b = wc*(1-pS)+wb*(pS);
        wbar_b = (x0*wbh_b+(N-x0)*wb)/N;
        p_b=(x0/N)*wbh_b/wbar_b;
        y_b = binopdf(x,N,p_b);
        if isnan(mean(y_b))
            y_b = binopdf(x,N,1.0);
        end
        yb = y_b.*phenodist(j); %weight transition probability vector by probability of having this realized pS
        tmatb = [tmatb; yb];
        
        
    end
    
    yya = sum(tmata); %weighted mean of all transition probabilities
    Ta = [Ta yya'];
    
    yyb = sum(tmatb); %weighted mean of all transition probabilities
    Tb = [Tb yyb'];
end

if emodel == 0 %if no env stochasticity, only 2 possible env regimes
    AB = Tb*Ta;
    BA = Ta*Tb;
    T = (AB+BA)./2;
elseif emodel == 1 %if env stochasticity, 4 possible env regimes
    AB = Tb*Ta;
    BA = Ta*Tb;
    AA = Ta*Ta;
    BB = Tb*Tb;
    T = (AB+BA+AA+BB)./4;
end
    

end
