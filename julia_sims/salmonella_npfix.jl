using Distributions
using TickTock

#warning, this code will take forever to run as is. parallelization is highly recommended

#environmental and phenotypic stochasticity are both turned "on"
envs = true
etype = "e1"
phenos = true
ptype = "p1"

filename = "salmonella_npfix.csv"
println(filename)

all_pA = [0.9, 0.931, 0.935, 0.94, 0.95, 0.975] #pA = probability of not taking antibitiotics
pSpec = 0.81

maxn = 9 #log10 maximum population size surveyed
lengthn = 10 #number of population sizes 
a=range(0, stop=maxn, length = lengthn)
allN=[convert(Int64,floor(10^i)) for i in a] #creates a vector of population sizes surveyed

tick()

function stoch_env_select(pA) #function that randomly generates the environment
    env = rand()
    return env <= pA
    #if returns true ==> in env. A; if returns false ==> in env. B
end

function det_env_select(generations, env_init) #function that deterministically generates the environment
    if env_init
        if isodd(generations)
            env = true
        else
            env = false
        end
    else
        if isodd(generations)
            env = false
        else
            env = true
        end
    end
end

function stoch_pheno(pSpec, BHcount) #function that randomly generates the phenotypic distribution of bet hedger offspring
    SpecPhe = rand(Binomial(BHcount, pSpec))
    #generates number of bet hedgers with specialist phenotype
    ConsPhe = BHcount-SpecPhe
    return [SpecPhe, ConsPhe]
    #returns a list with the # of bet hedgers with specialist and conservative phenotypes respectively
end

function det_pheno(pSpec, BHcount) #function that deterministically calculates the phenotypic distribution of bet hedger offspring
    SpecPhe = BHcount*pSpec
    #generates number of bet hedgers with specialist phenotype
    ConsPhe = BHcount-SpecPhe
    return [SpecPhe, ConsPhe]
    #returns a list with the # of bet hedgers with specialist and conservative phenotypes respectively
end

function average_fitness(PopNum, wBH, wWT)
    return (PopNum[1]*wBH+PopNum[2]*wWT)/(sum(PopNum))
end

function reproduction(env, PopNum, PopSize, spec, cons, phenotype_counts)
    #simulates random wright fisher reproduction in one generation
    wBH = Float64 #bet hedger fitness
    wWT = Float64 #wild type specialist fitness
    if env 
        #if in env. A
        wBH = phenotype_counts[1]/PopNum[1]*spec[1]+phenotype_counts[2]/PopNum[1]*cons[1]
        wWT = spec[1] #WT always has fast-growing phenotype
    else
        #else in env. B
        wBH = phenotype_counts[1]/PopNum[1]*spec[2]+phenotype_counts[2]/PopNum[1]*cons[2]
        wWT = spec[2]
    end
    wBar = average_fitness(PopNum, wBH, wWT)
    p = PopNum[1]/PopSize * wBH/wBar
    PopNum[1] = rand(Binomial(PopSize, p))
    #number of bet hedger offspring is pulled randomly from a binomial distribution whose
    #expectation is based on the relative fitness of the bet hedger
    PopNum[2] = PopSize-PopNum[1]
    return PopNum
end

function generation(PopNum, spec, cons, pA, env_init, pSpec, generations, envs, phenos)
    #determine the environment for this generation
    if envs
        env = stoch_env_select(pA)
    else
        env = det_env_select(generations, env_init)
    end
    #determine the realized phenotypic makeup of the BH
    if phenos
        phenotype_counts = stoch_pheno(pSpec, PopNum[1]) 
    else
        phenotype_counts = det_pheno(pSpec, PopNum[1])
    end
    #simulate random wright fisher reproduction based on fitness values given by environment and phenotypic distribution
    PopNum = reproduction(env, PopNum, sum(PopNum), spec, cons, phenotype_counts)
    return PopNum
end

function simulate(PopSize::Int64, pA::Float64, pSpec::Float64, envs::Bool, phenos::Bool)
    InitNum = 1 #number of bh to start
    PopNum = [InitNum, PopSize-InitNum]
    #PopNum = [# of BH, # of WT]
    spec = [1.8, 0.003] #WT fast growing / antibiotic sensitive phenotype in Env. A (no antibiotics) and Env. B (yes antibiotics) respectively
    cons = [1, 0.04] #persister slow growing / antibiotic tolerant phenotype in Env. A (no antibiotics) and Env. B (yes antibiotics) respectively
    generations::Int64 = 1
    env_init = convert(Bool, rand(Binomial(1,pA)))
    while 0<PopNum[1]<PopSize
        PopNum = generation(PopNum, spec, cons, pA, env_init, pSpec, generations, envs, phenos)
        generations+=1
    end
    return PopNum[1]==PopSize
    #returns True (1) if BH reaches fixations, False (0) if BH is lost
end

colnames = ["N", "NPfix", "pA"]
out = open(filename, "w") #creates a new output file 
write(out, join(colnames, ","), "\n") #populates output file with column names
close(out)

for pA in all_pA
    c = 0 #counts number of replicates that reach fixation
    for N in allN
        if N < 10^5
            reps = 1000*10^5
        else
            reps = 100*N
        end
        for run = 1:reps
            c += simulate(N, pA, pSpec, envs, phenos) #c increases by 1 for each rep that reaches fixation
        end
        npf = ((c/reps)*N)
        output = [N, npf, pA]
        out = open(filename, "a") #adds the NPfix value to the output file
        write(out, join(output, ","), "\n")
        close(out)

        println(npf)
    end
end

tock()
