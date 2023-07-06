function T = e0_meanmat(wa, wb, s, pS, N)
    g = 1;
    gmf = geomean([wa, wb]) + s;
    wc = (5*pS - 2*(4*gmf^2 + (9*pS^2)/4)^(1/2))/(4*pS - 4); 
    
    wbha = wc*(1-pS)+wa*(pS);
    Ta = construct_markov(N, wbha, wa);
    wbhb = wc*(1-pS)+wb*(pS);
    Tb = construct_markov(N, wbhb, wb);
    
    AB = Tb^g*Ta^g;
    BA = Ta^g*Tb^g;
    
    T = (AB+BA)./2;
end