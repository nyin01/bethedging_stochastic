using Distributions
using TickTock

if ARGS[1] == "1"
    envs = true
    etype = "e1"
else
    envs = false
    etype = "e0"
end

if ARGS[2] == "1"
    phenos = true
    ptype = "p1"
else
    phenos = false
    ptype = "p0"
end

#pA = 0.5 #probability of being in env. A
#all_pA = [0.09, 0.118, 0.15, 0.2, 0.396, 0.425]
#pSpec = parse(Float64, ARGS[3]) #proportion of BH offspring with specialist phenotype
pSpec = 0.53

wc = 1.0

filename = "821papaver_fall_wt_0.836$etype.$ptype.csv"
println(filename)

#maxn = parse(Int64, ARGS[4])
maxn = 5
lengthn = 10
a=range(0, stop=maxn, length = lengthn)
allN=[convert(Int64,floor(10^i)) for i in a] #creates a vector of 10 population
#sizes that are evenly spaced on a log scale


reps = floor(1000*10^maxn) #number of replicates

tick()

function stoch_env_select(pA) #function that randomly generates the environment
    env = rand()
    return env <= pA
    #if returns true ==> in env. A; if returns false ==> in env. B
end

function det_env_select(generations, env_init)
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

function stoch_pheno(pSpec, BHcount)
    SpecPhe = rand(Binomial(BHcount, pSpec))
    #generates number of bet hedgers with specialist phenotype
    ConsPhe = BHcount-SpecPhe
    return [SpecPhe, ConsPhe]
    #returns a list with the # of bet hedgers with specialist and conservative phenotypes respectively
end

function det_pheno(pSpec, BHcount)
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
    if env #collapse if else statement
        #if in env. A
        wBH = phenotype_counts[1]/PopNum[1]*spec[1]+phenotype_counts[2]/PopNum[1]*cons[1]
        wWT = spec[1]
        #wWT = cons[1]
    else
        #else in env. B
        wBH = phenotype_counts[1]/PopNum[1]*spec[2]+phenotype_counts[2]/PopNum[1]*cons[2]
        wWT = spec[2]
        #wWT = cons[2]
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
    if envs
        env = stoch_env_select(pA) #determines the environment for this generation
    else
        env = det_env_select(generations, env_init)
    end
    if phenos
        phenotype_counts = stoch_pheno(pSpec, PopNum[1]) #determines the realized phenotypic makeup of the BH
    else
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
    spec = [18, 0.836] #specialist phenotype in Env. A and Env. B respectively
    cons = [wc, wc] #conservative phenotype in Env. A and Env. B respectively
    generations::Int64 = 1
    env_init = convert(Bool, rand(Binomial(1,pA)))
    while 0<PopNum[1]<PopSize
        PopNum = generation(PopNum, spec, cons, pA, env_init, pSpec, generations, envs, phenos)
        generations+=1
    end
    return PopNum[1]==PopSize
    #returns True (1) if BH win, False (0) otherwise
end

println(allN)

colnames = ["N", "NPfix", "pA"]
out = open(filename, "w") #creates a new output file whose filename includes parameters
write(out, join(colnames, ","), "\n") #populates output file with vector of 10 population sizes
close(out)

for pA in all_pA
    Npfix = Float64[] #creates an empty vector where normalized pfix values will be added
    for N in allN #repeats this process at each population size
        c = 0 #counts number of replicates that reach fixation
        for run = 1:reps
            c += simulate(N, pA, pSpec, wc, envs, phenos) #c increases by 1 for each rep that reaches fixation
        end
        push!(Npfix,((c/reps)*N))
        output = [N, ((c/reps)*N), pA]
        out = open(filename, "a") #adds the NPfix vector to the output file
        write(out, join(output, ","), "\n")
        close(out)
    end
    println(Npfix)
end

tock()
