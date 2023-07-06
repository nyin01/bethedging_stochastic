Nsizes = [1, 2, 4, 7, 11, 19, 31, 51, 84, 138, 227, 372];

g = 1;

wa = 2;       %wt in env. a
wb = 0.5;     %wt in env. b
pS = 0.5;
MaxGen = 1000;
s = 0;

outm = [];

% e1p1 model
for j = 1:length(Nsizes)
    N = Nsizes(j);
    mat = p1_meanmat(wa, wb, s, pS, N, 1);
    npf = markov_multiplication(mat, MaxGen);
    outm = [outm; N, npf, "e1p1"];
end

% e1p0 model
for j = 1:length(Nsizes)
    N = Nsizes(j);
    mat = e1_meanmat(wa,wb,s,pS,N);
    npf = markov_multiplication(mat, MaxGen);
    outm = [outm; N, npf, "e1p0"];
end

% e0p1 model
for j = 1:length(Nsizes)
    N = Nsizes(j);
    mat = p1_meanmat(wa, wb, s, pS, N, 0);
    npf = markov_multiplication(mat, MaxGen);
    outm = [outm; N, npf, "e0p1"];
end

% e0p0 model
for j = 1:length(Nsizes)
    N = Nsizes(j);
    mat = e0_meanmat(wa, wb, s, pS, N);
    npf = markov_multiplication(mat, MaxGen);
    outm = [outm; N, npf, "e0p0"];
end

    
out = array2table(outm); 
out.Properties.VariableNames(1:3) = {'N', 'npf', 'stoch'}; 
writetable(out, "markovnpfix_stochtype.csv");
