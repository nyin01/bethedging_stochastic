#code to calculate normalized probability of fixation curves

using Distributions
# using TickTock

# ============ ARGS ============
#1) environmental stochasticity (on or off)
#2) germination stochasticity (on or off)
#3) soft/hard selection
#4) probability of germination
#5) seed bank length
fits = parse(Bool, ARGS[1])
ftype = fits ? "fs" : "fd"
germs = parse(Bool, ARGS[2])
gtype = germs ? "gs" : "gd"
hards = parse(Bool, ARGS[3])
stype = hards ? "hs" : "ss"
pGermBH = parse(Float64, ARGS[4])
pGermWT = parse(Float64, ARGS[5])
bankLength = parse(Int64, ARGS[6]) #generations seeds are allowed to stay dormant in the seed bank before being discarded
# fitness dist params
fmean = parse(Float64, ARGS[7])
fvar = parse(Float64, ARGS[8])
# site label
site = ARGS[9]

filename = "delayedgerm_$ftype:$gtype:$stype:$pGermBH:$pGermWT:$bankLength:$fmean:$fvar:$site.csv"
save_dir = "data/$site" * "_dist/"
save_path = save_dir * filename

# ============ META PARAMS ============
maxn = 3 #10^maxn is the largest population size
lengthn = 10 #number of population sizes to repeat
a=range(0, stop=maxn, length=lengthn)
allN=[convert(Int64,floor(10^i)) for i in a] #creates a vector of 10 population sizes that are evenly spaced on a log scale
reps = floor(100*10^maxn) #number of replicates

# ============ FUNCTIONS ============
function stoch_fitness_select(mean, var)
    #function that randomly generates the fitness value
    w = rand(Normal(mean, sqrt(var)))
    return w < 0 ? 0 : w # fitness cannot be negative
end

function det_fitness_select(mean, var, generations)
    #function that deterministically generates the fitness value
    w = 0.0
    if isodd(generations)
        w = mean + sqrt(var)
    else
        w = mean - sqrt(var)
    end
    return w < 0 ? 0 : w # fitness cannot be negative
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

function hard_reproduction(PopNum, K, w)
    # Check if the sum of the population is zero
    if sum(PopNum) == 0
        return [0, 0]
    end

    # Calculate new population numbers using the Poisson distribution
    lambdaBH = w * PopNum[1]
    lambdaWT = w * PopNum[2]

    PopNum = [rand(Poisson(lambdaBH)), rand(Poisson(lambdaWT))]

    # Check if the total population exceeds the carrying capacity
    if sum(PopNum) > K
        # Calculate proportion and ensure it is within [0,1]
        p = PopNum[1] / sum(PopNum)
        # Assign new population numbers based on binomial distribution
        PopNum[1] = rand(Binomial(K, p))
        PopNum[2] = K - PopNum[1]
    end

    return PopNum
end


function soft_reproduction(PopNum, PopSize)
    # Check if the sum of the population is zero
    if sum(PopNum) == 0
        return [0, 0]
    end

    # Calculate the proportion of bet hedgers
    p = PopNum[1] / sum(PopNum)

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

function generation(PopNum, K, generations, BHbank, WTbank, pGermBH, pGermWT, germs, fits, fmean, fvar)
    #generation begins by determining the fitness
    w = fits ? stoch_fitness_select(fmean, fvar) : det_fitness_select(fmean, fvar, generations)
    #generate phenotypic distribution amongst diversified BH offspring
    PopNum, BHbank, WTbank = germination(PopNum, BHbank, WTbank, pGermBH, pGermWT, germs)
    # hard selection: use sum of curent population size; soft selection: use fixed carrying capacity
    PopNum = hards ? hard_reproduction(PopNum, K, w) : soft_reproduction(PopNum, K)
    #simulates random wright fisher reproduction
    return PopNum, BHbank, WTbank
end

function simulate(PopSize::Int64, bankLength::Int64, fits::Bool, germs::Bool, fmean::Float64, fvar::Float64)
    InitNum = 1 #number of bh active to start
    BHbank = zeros(Int64, bankLength) #creates a vector of zeros of length bankLength
    WTbank = zeros(Int64, bankLength) #creates a vector of zeros of length bankLength
    PopNum = [InitNum, PopSize-InitNum] #PopNum = [# of BH active, # of WT active]
    generations::Int64 = 1
    # continue the simulation when the BH is not lost or the BHbank is not empty
    while (0<PopNum[1]<PopSize || sum(BHbank)>0) && generations<PopSize*100
        #while the BH has not been lost (0%) or reached fixation (100%)
        PopNum, BHbank, WTbank = generation(PopNum, PopSize, generations, BHbank, WTbank, pGermBH, pGermWT, germs, fits, fmean, fvar)
        generations+=1
    end
    return PopNum[1]>=PopSize
    #returns True (1) if BH reaches fixation, False (0) if the BH is lost
end

# ============ SIMULATION ============
#begin
# tick()
println("ftype: ", ftype)
println("gtype: ", gtype)
println("stype: ", stype)
println("pGermBH: ", pGermBH)
println("pGermWT: ", pGermWT)
println("bankLength: ", bankLength)
println("fmean: ", fmean)
println("fvar: ", fvar)
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
        c += simulate(N, bankLength, fits, germs, fmean, fvar) #c increases by 1 for each rep that reaches fixation
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
