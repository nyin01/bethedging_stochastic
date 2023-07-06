function T = e1_meanmat(wa, wb, s, pS, N)
    %effective generation matrix for E1P0 model
    g = 1;
    gmf = geomean([wa, wb]) + s;
    wc = (5*pS - 2*(4*gmf^2 + (9*pS^2)/4)^(1/2))/(4*pS - 4); 
    
    wbha = wc*(1-pS)+wa*(pS);
    Ta = construct_markov(N, wbha, wa); %env a matrix if pSpec is exactly equal to the expectation
    wbhb = wc*(1-pS)+wb*(pS);
    Tb = construct_markov(N, wbhb, wb); %env b matrix if pSpec is exactly equal to the expectation

    %four possible environmental regimes
    AB = Tb^g*Ta^g;
    BA = Ta^g*Tb^g;
    AA = Ta^g*Ta^g;
    BB = Tb^g*Tb^g;
    
    T = (AB+BA+AA+BB)./4;
end
