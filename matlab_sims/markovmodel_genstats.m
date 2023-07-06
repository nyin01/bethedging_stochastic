%calculate stats (ploss and expected offspring) from markov model for the first generation and effective
%generation

outm = []; %output matrix

Nsizes = [100]; %population size
wa = 2;       %wt fitness in env. a
wb = 0.5;     %wt fitness in env. b
pSpec = 0.5; %probability of adopting specialist phenotype for BH
s = 0; %delta gmf
gmf = 1 + s;
wc = (5*pSpec - 2*(4*gmf^2 + (9*pSpec^2)/4)^(1/2))/(4*pSpec - 4); %calculates wC as function of deltaGMF and pSpec
wbha = wa*pSpec+wc*(1-pSpec); %expected fitness of BH in A
wbhb = wb*pSpec+wc*(1-pSpec); %expected fitness of BH in B

%step 1: first generation
for N = Nsizes
    x = [0:N];
    %in the first generation, there is only one individual. therefore there
    %is a finite number of possibilities (2) for the phenotypic
    %distribution
    
    Ta = construct_markov(N, wbha, wa); %matrix for env a if BH has exactly expected pSpec
    apL = Ta(1,2); %probability of loss --> probability x(t+1) = 0 given x(t) = 1
    aW = sum(Ta(:,2).*x'); %expected offspring --> given x(t) = 1, what is the mean # of individuals in x(t+1)
    
    Tb = construct_markov(N, wbhb, wb); %matrix for env b if BH has exactly expected pSpec
    bpL = Tb(1,2); %probability of loss
    bW = sum(Tb(:,2).*x'); %expected offspring
    
    Tac = construct_markov(N, wc, wa); %matrix for env a if BH is all conservative phenotype
    acpL = Tac(1,2); %probability of loss
    acW = sum(Tac(:,2).*x'); %expected offspring
    
    Tbc = construct_markov(N, wc, wb); %matrix for env b if BH is all conservative phenotype
    bcpL = Tbc(1,2); %probability of loss
    bcW = sum(Tbc(:,2).*x'); %expected offspring
    
    Tn = construct_markov(N, gmf, gmf); %matrix for if the BH is neutral
    %i.e. if the BH adopts the specialist phenotype in either env A or B
    npL = Tn(1,2); %probability of loss
    nW = sum(Tn(:,2).*x'); %expected offspring
    
    p1_pL = mean([acpL, bcpL, npL, npL]); %mean pLoss for the P1 model
    p1_W = mean([acW, bcW, nW, nW]); %mean expected offspring for the P1 model
    p0_pL = mean([apL, bpL]); %mean pLoss for the P0 model
    p0_W = mean([aW, bW]); %mean expected offspring for the P0 model
    
    %because env stochasticity does not matter in first generation, only
    %pheno stoch type is used
    outm = [outm; N, "E1P1", "g", p1_pL, p1_W]; 
    outm = [outm; N, "E1P0", "g", p0_pL, p0_W];
    outm = [outm; N, "E0P1", "g", p1_pL, p1_W];
    outm = [outm; N, "E0P0", "g", p0_pL, p0_W];
    
end


%step 2 = effective generation
for N = Nsizes
    x = [0:N];
    e1p1 = p1_meanmat(wa, wb, s, pSpec, N, 1); %constructs matrix for E1P1
    e1p1_L = e1p1(1,2); %probability of loss 
    e1p1_W = sum(e1p1(:,2).*x'); %expected offspring
    
    e1p0 = e1_meanmat(wa, wb, s, pSpec, N); %constructs matrix for E1P0
    e1p0_L = e1p0(1,2);
    e1p0_W = sum(e1p0(:,2).*x');
    
    e0p1 = p1_meanmat(wa, wb, s, pSpec, N, 0); %constructs matrix for E0P1
    e0p1_L = e0p1(1,2);
    e0p1_W = sum(e0p1(:,2).*x');
    
    e0p0 = e0_meanmat(wa, wb, s, pSpec, N); %constructs matrix for E0P0
    e0p0_L = e0p0(1,2);
    e0p0_W = sum(e0p0(:,2).*x');
    
    %pushes stats to output matrix
    outm = [outm; N, "E1P1", "eg", e1p1_L, e1p1_W];
    outm = [outm; N, "E1P0", "eg", e1p0_L, e1p0_W];
    outm = [outm; N, "E0P1", "eg", e0p1_L, e0p1_W];
    outm = [outm; N, "E0P0", "eg", e0p0_L, e0p0_W];
end

%writes to table
out = array2table(outm); 
out.Properties.VariableNames(1:5) = {'N', 'stoch', 'gentype', 'pLoss', 'w'}; 
