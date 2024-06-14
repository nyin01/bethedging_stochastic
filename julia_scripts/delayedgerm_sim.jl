#code to calculate normalized probability of fixation curves

using Distributions
# using TickTock

# ============ ARGS ============
#1) environmental stochasticity (on or off)
#2) germination stochasticity (on or off)
#3) soft/hard selection
#4) probability of germination
#5) seed bank length
envs = parse(Bool, ARGS[1])
etype = envs ? "es" : "ed"
germs = parse(Bool, ARGS[2])
gtype = germs ? "gs" : "gd"
hards = parse(Bool, ARGS[3])
stype = hards ? "hs" : "ss"
pA = parse(Float64, ARGS[4]) #probability of being in env. A
pGermBH = parse(Float64, ARGS[5])
pGermWT = parse(Float64, ARGS[6])
bankLength = parse(Int64, ARGS[7]) #generations seeds are allowed to stay dormant in the seed bank before being discarded
wA = parse(Float64, ARGS[8]) #fitness of the bet hedger in env. A
wB = parse(Float64, ARGS[9]) #fitness of the bet hedger in env. B
site = ARGS[10]
w = [wA, wB] #specialist phenotype in Env. A and Env. B respectively. Given pA = 0.5 yields a GMF of 1
filename = "delayedgerm_$etype:$gtype:$stype:$pA:$pGermBH:$pGermWT:$bankLength:$wA:$wB:$site.csv"
save_dir = "data/June13/$site" * "_ab/"
save_path = save_dir * filename

# ============ META PARAMS ============
maxn = 4 #10^maxn is the largest population size
lengthn = 10 #number of population sizes to repeat
a=range(0, stop=maxn, length=lengthn)
allN=[convert(Int64,floor(10^i)) for i in a] #creates a vector of 10 population sizes that are evenly spaced on a log scale
reps = floor(100*10^maxn) #number of replicates

# ============ FUNCTIONS ============
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

function stoch_germ_select(seeds::Int64, pGerm::Float64)
    #function that randomly germinates seeds
    #returns [germinated, dormant]
    germinated = rand(Binomial(seeds, pGerm))
    return germinated, seeds-germinated
end

function det_germ_select(seeds::Int64, pGerm::Float64)
    #function that deterministically germinates seeds
    #returns [germinated, dormant]
    germinated = round(seeds*pGerm)
    return germinated, seeds-germinated
end

# function hard_reproduction(env, PopNum, K, w)
#     if sum(PopNum) == 0
#         return [0, 0]
#     end
#     #otherwise
#     #simulates random wright fisher reproduction in one generation
#     wBH = Float64 #bet hedger fitness
#     wWT = Float64 #wild type fitness
#     #in this model, fitness is the same for both phenotypes in the same environment
#     if env
#         #if in env. A
#         wBH = w[1]
#         wWT = w[1]
#     else
#         #else in env. B
#         wBH = w[2]
#         wWT = w[2]
#     end
#     PopNum = [rand(Poisson(wBH*PopNum[1])), rand(Poisson(wWT*PopNum[2]))]
#     if sum(PopNum) > K
#         p = PopNum[1]/sum(PopNum)
#         PopNum[1] = rand(Binomial(K, p))
#         #number of bet hedger offspring is pulled randomly from a binomial distribution whose
#         #expectation is based on the relative fitness of the bet hedger
#         PopNum[2] = K-PopNum[1]
#     end
#     return PopNum
# end

function hard_reproduction(env, PopNum, K, w)
    # Check if the sum of the population is zero
    if sum(PopNum) == 0
        return [0, 0]
    end

    # Ensure PopNum values are non-negative integers
    if any(x -> x < 0 || x != floor(x), PopNum)
        throw(ArgumentError("Population numbers must be non-negative integers."))
    end

    # Ensure K is a positive integer
    if K <= 0 || K != floor(K)
        throw(ArgumentError("Carrying capacity K must be a positive integer."))
    end

    # Initialize fitness values
    wBH = 0.0 # bet hedger fitness
    wWT = 0.0 # wild type fitness

    # Assign fitness values based on the environment
    if env
        wBH = w[1]
        wWT = w[1]
    else
        wBH = w[2]
        wWT = w[2]
    end

    # Ensure fitness values lead to valid Poisson parameters
    if wBH < 0 || wWT < 0
        throw(ArgumentError("Fitness values must be non-negative."))
    end

    # Calculate new population numbers using the Poisson distribution
    lambdaBH = wBH * PopNum[1]
    lambdaWT = wWT * PopNum[2]

    if lambdaBH < 0 || lambdaWT < 0
        throw(ArgumentError("Calculated Poisson parameters must be non-negative."))
    end

    PopNum = [
        rand(Poisson(lambdaBH)),
        rand(Poisson(lambdaWT))
    ]

    # Check if the total population exceeds the carrying capacity
    if sum(PopNum) > K
        # Handle division by zero
        total_population = sum(PopNum)
        if total_population == 0
            return [0, 0]
        end

        # Calculate proportion and ensure it is within [0,1]
        p = PopNum[1] / total_population
        p = clamp(p, 0.0, 1.0)

        # Assign new population numbers based on binomial distribution
        PopNum[1] = rand(Binomial(K, p))
        PopNum[2] = K - PopNum[1]
    end

    return PopNum
end


# function reproduction(env, PopNum, PopSize, w)
#     #if no active population, return [0, 0]
#     #this function uses rel fitness
#     if sum(PopNum) == 0
#         return [0, 0]
#     end
#     #otherwise
#     #simulates random wright fisher reproduction in one generation
#     wBH = Float64 #bet hedger fitness
#     wWT = Float64 #wild type fitness
#     #in this model, fitness is the same for both phenotypes in the same environment
#     if env
#         #if in env. A
#         wBH = w[1]
#         wWT = w[1]
#     else
#         #else in env. B
#         wBH = w[2]
#         wWT = w[2]
#     end
#     p = PopNum[1]/PopSize
#     PopNum[1] = rand(Binomial(PopSize, p))
#     #number of bet hedger offspring is pulled randomly from a binomial distribution whose
#     #expectation is based on the relative fitness of the bet hedger
#     PopNum[2] = PopSize-PopNum[1]
#     #returns vector of [# of BH offspring, # of WT offspring]
#     return PopNum
# end

function reproduction(env, PopNum, PopSize, w)
    # Check if the sum of the population is zero
    if sum(PopNum) == 0
        return [0, 0]
    end

    # Ensure PopNum values are non-negative integers
    if any(x -> x < 0 || x != floor(x), PopNum)
        throw(ArgumentError("Population numbers must be non-negative integers."))
    end

    # Ensure PopSize is a positive integer
    if PopSize <= 0 || PopSize != floor(PopSize)
        throw(ArgumentError("Population size must be a positive integer."))
    end

    # Initialize fitness values
    wBH = 0.0 # bet hedger fitness
    wWT = 0.0 # wild type fitness

    # Assign fitness values based on the environment
    if env
        # if in env. A
        wBH = w[1]
        wWT = w[1]
    else
        # else in env. B
        wBH = w[2]
        wWT = w[2]
    end

    # Calculate the proportion of bet hedgers
    p = PopNum[1] / sum(PopNum)

    if p < 0 || p > 1
        # println("p: ", p)
        throw(ArgumentError("Proportion of bet hedgers must be within [0,1]."))
    end

    # Assign new population numbers based on binomial distribution
    PopNum[1] = rand(Binomial(PopSize, p))
    PopNum[2] = PopSize - PopNum[1]

    # Returns vector of [# of BH offspring, # of WT offspring]
    return PopNum
end


function germination(PopNum, BHbank, WTbank, pGermBH, pGermWT, germs)
    #choose the germination function
    fGerm = germs ? stoch_germ_select : det_germ_select
    #all candidate seeds, including the active population and the seed bank, length=bankLength+1
    BHcandidates = [PopNum[1]; BHbank]
    WTcandidates = [PopNum[2]; WTbank]
    #map germination function to the vector
    BHtuples = fGerm.(BHcandidates, pGermBH) #BH uses a germination rate
    WTtuples = fGerm.(WTcandidates, pGermWT) #WT has 100% germination (using the same function for consistency)
    #separate germinated and dormant seeds
    BHgermed = [t[1] for t in BHtuples]
    BHdormant = [t[2] for t in BHtuples]
    WTgermed = [t[1] for t in WTtuples]
    WTdormant = [t[2] for t in WTtuples]
    #update PopNum with germinated seeds
    PopNum[1] = sum(BHgermed)
    PopNum[2] = sum(WTgermed)
    #deposit dormant seeds back into the seed bank, shifting one generation and discarding the oldest generation
    BHbank = BHdormant[1:end-1]
    WTbank = WTdormant[1:end-1]
    #return updated PopNum, BHbank, WTbank
    return PopNum, BHbank, WTbank
end

function generation(PopNum, K, w, pA, env_init, generations, envs, BHbank, WTbank, pGermBH, pGermWT, germs)
    #generation begins by determining the environment
    env = envs ? stoch_env_select(pA) : det_env_select(generations, env_init)
    #generate phenotypic distribution amongst diversified BH offspring
    PopNum, BHbank, WTbank = germination(PopNum, BHbank, WTbank, pGermBH, pGermWT, germs)
    # hard selection: use sum of curent population size; soft selection: use fixed carrying capacity
    PopNum = hards ? hard_reproduction(env, PopNum, K, w) : reproduction(env, PopNum, K, w)
    #simulates random wright fisher reproduction
    return PopNum, BHbank, WTbank
end

function simulate(PopSize::Int64, pA::Float64, bankLength::Int64, envs::Bool, germs::Bool, w::Vector{Float64})
    InitNum = 1 #number of bh active to start
    BHbank = zeros(Int64, bankLength) #creates a vector of zeros of length bankLength
    WTbank = zeros(Int64, bankLength) #creates a vector of zeros of length bankLength
    PopNum = [InitNum, PopSize-InitNum] #PopNum = [# of BH active, # of WT active]
    generations::Int64 = 1
    env_init = convert(Bool, rand(Binomial(1,pA))) #ensures first generation is random for E0 model
    # continue the simulation when the BH is not lost or the BHbank is not empty
    while (0<PopNum[1]<PopSize || sum(BHbank)>0) && generations<PopSize*100
        #while the BH has not been lost (0%) or reached fixation (100%)
        PopNum, BHbank, WTbank = generation(PopNum, PopSize, w, pA, env_init, generations, envs, BHbank, WTbank, pGermBH, pGermWT, germs)
        generations+=1
    end
    return PopNum[1]>=PopSize
    #returns True (1) if BH reaches fixation, False (0) if the BH is lost
end

# ============ SIMULATION ============
#begin
# tick()
println("etype: ", etype)
println("gtype: ", gtype)
println("stype: ", stype)
println("pA: ", pA)
println("pGermBH: ", pGermBH)
println("pGermWT: ", pGermWT)
println("bankLength: ", bankLength)
println("wA: ", wA)
println("wB: ", wB)
println("site: ", site)
println("filename: ", filename)
colnames = ["N", "NPfix"]
global out = open(save_path, "w") #creates a new output file whose filename includes parameters
write(out, join(colnames, ","), "\n")
close(out)
Npfix = Float64[] #creates an empty vector where normalized pfix values will be added
for N in allN #repeats this process at each population size
    c = 0 #counts number of replicates that reach fixation
    for run = 1:reps
        c += simulate(N, pA, bankLength, envs, germs, w) #c increases by 1 for each rep that reaches fixation
    end
    push!(Npfix,((c/reps)*N))
    output = [N, ((c/reps)*N)]
    global out = open(save_path, "a") #adds the NPfix vector to the output file
    write(out, join(output, ","), "\n")
    close(out)
end
#end
# tock()
println(Npfix)
