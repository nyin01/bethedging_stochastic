#code to calculate normalized probability of fixation curves

using Distributions
using TickTock

#takes in 4 arguments: 1) environmental stochasticity (on or off); 2) phenotypic stochasticity (on or off); 3) pSpec of the invading bet hedger; 4) Delta GMF
if ARGS[1] == "1"
    envs = true
    etype = "es"
else
    envs = false
    etype = "ed"
end

if ARGS[2] == "1"
    phenos = true
    ptype = "ps"
else
    phenos = false
    ptype = "pd"
end

pA = 0.5 #probability of being in env. A
pSpec = parse(Float64, ARGS[3]) #proportion of BH offspring with specialist phenotype
#if 0, conservative bet hedger. if 1, WT specialist. if float, diversified bet hedger
#gmf = parse(Float64, ARGS[3])
s = parse(Float64, ARGS[4]) #Delta GMF, expected difference between fitness of the bet hedger and the wild type
gmf = 1+s
wc = round((5*pSpec - 2*(4*gmf^2 + (9*pSpec^2)/4)^(1/2))/(4*pSpec - 4), digits = 3)#calculates WC for the bet hedger given pSpec and GMF

filename = "bhfix_$gmf.$pSpec.$etype.$ptype.csv" #output file
println(filename)

maxn = 5 #10^maxn is the largest population size
lengthn = 20 #number of population sizes to repeat
a=range(0, stop=maxn, length = lengthn)
allN=[convert(Int64,floor(10^i)) for i in a] #creates a vector of 10 population sizes that are evenly spaced on a log scale

reps = floor(1000*10^maxn) #number of replicates

tick()

function stoch_env_select(pA)
    #function that randomly generates the environment, E1 model
    env = rand()
    return env <= pA
    #if returns true ==> in env. A; if returns false ==> in env. B
end

function det_env_select(generations, env_init)
    #function that deterministically switches the environment every generation, E0 model
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
    #if returns true ==> in env. A; if returns false ==> in env. B
end

function stoch_pheno(pSpec, BHcount)
    #function that stochastically generates phenotypic distribution of the BH, P1 model
    SpecPhe = rand(Binomial(BHcount, pSpec)) #randomly pulls number of bet hedgers with specialist phenotype
    ConsPhe = BHcount-SpecPhe #remainder adopt conservative phenotype
    return [SpecPhe, ConsPhe]
    #returns a list with the # of bet hedgers with specialist and conservative phenotypes respectively
end

function det_pheno(pSpec, BHcount)
    #function that deterministically generates phenotypic distribution of the BH, P0 model
    SpecPhe = BHcount*pSpec #number of bet hedgers with specialist phenotype is exactly the expectation
    ConsPhe = BHcount-SpecPhe #remainder adopt conservative phenotype
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
    if env #collapse if else statement
        #if in env. A
        wBH = phenotype_counts[1]/PopNum[1]*spec[1]+phenotype_counts[2]/PopNum[1]*cons[1]
        wWT = spec[1]
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
    #returns vector of [# of BH offspring, # of WT offspring]
    return PopNum
end

function generation(PopNum, spec, cons, pA, env_init, pSpec, generations, envs, phenos)
    #generation begins by determining the environment
    if envs #if E1 model
        env = stoch_env_select(pA) #determines the environment for this generation
    else #elseif E0 model
        env = det_env_select(generations, env_init)
    end
    #generate phenotypic distribution amongst diversified BH offspring
    if phenos #if P1 model
        phenotype_counts = stoch_pheno(pSpec, PopNum[1]) #determines the realized phenotypic makeup of the BH
    else #elseif P0 model
        phenotype_counts = det_pheno(pSpec, PopNum[1])
    end
    PopNum = reproduction(env, PopNum, sum(PopNum), spec, cons, phenotype_counts)
    #simulates random wright fisher reproduction
    return PopNum
end

function simulate(PopSize::Int64, pA::Float64, pSpec::Float64, wc::Float64, envs::Bool, phenos::Bool)
    InitNum = 1 #number of bh to start
    PopNum = [InitNum, PopSize-InitNum]
    #PopNum = [# of BH, # of WT]
    spec = [2, 0.5] #specialist phenotype in Env. A and Env. B respectively. Given pA = 0.5 yields a GMF of 1
    cons = [wc, wc] #conservative phenotype in Env. A and Env. B respectively
    generations::Int64 = 1
    env_init = convert(Bool, rand(Binomial(1,pA))) #ensures first generation is random for E0 model
    while 0<PopNum[1]<PopSize
        #while the BH has not been lost (0%) or reached fixation (100%)
        PopNum = generation(PopNum, spec, cons, pA, env_init, pSpec, generations, envs, phenos)
        generations+=1
    end
    return PopNum[1]==PopSize
    #returns True (1) if BH reaches fixation, False (0) if the BH is lost
end

println(allN)

colnames = ["N", "NPfix"]
out = open(filename, "w") #creates a new output file whose filename includes parameters
write(out, join(colnames, ","), "\n")
close(out)

Npfix = Float64[] #creates an empty vector where normalized pfix values will be added
for N in allN #repeats this process at each population size
    c = 0 #counts number of replicates that reach fixation
    for run = 1:reps
        c += simulate(N, pA, pSpec, wc, envs, phenos) #c increases by 1 for each rep that reaches fixation
    end
    push!(Npfix,((c/reps)*N))
    output = [N, ((c/reps)*N)]
    out = open(filename, "a") #adds the NPfix vector to the output file
    write(out, join(output, ","), "\n")
    close(out)
end
println(Npfix)
#end

tock()
