#extend npfix simulations across parameter space, automate classification of sign of selection

using Distributions
using Interpolations
using Dates
using TickTock
using Base.Threads

pA = 0.5 #probability of being in env. A
allpSpec = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] #proportion of BH offspring with specialist phenotype
allS = [0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075, 1.1, 1.125] #all values of GMF of the BH

maxn = 6
a=range(0, stop=maxn, length = 10)
allN=[convert(Int64,floor(10^i)) for i in a]

envs = true #environmental stochasticity on, E1
etype = "es"

phenos = true #phenotypic stochasticity on, P1
ptype = "ps"

reps = 500*10^maxn

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

function signflip1(vec) #for monotonic classifier, counts the number of times the derivative of the npfix curve changes sign
    sflip = 0
    for j in collect(2:length(vec))
        if vec[j]>0 && vec[j-1]<0
            sflip+=1
        elseif vec[j]<0 && vec[j-1]>0
            sflip+=1
        else
            sflip+=0
        end
    end
    return(sflip)
end

function signflip(vec) #for location classifier, counts the number of times the sign of selection changes
    sflip = 0
    for j in collect(2:length(vec))
        if vec[j]>=1 && vec[j-1]<1
            sflip+=1
        elseif vec[j]<1 && vec[j-1]>=1
            sflip+=1
        else
            sflip+=0
        end
    end
    return(sflip)
end

function monotonic(npfix)
    #classifies sign based on whether the slope of the npfix curve is monotonic
    difference = []
    for i in collect(2:length(npfix))
        push!(difference, npfix[i]-npfix[i-1])
    end
        if count(i->i==0, difference)==length(difference) #slope is always be 0
            Sign = 0 #neutral
        elseif count(i->i>=0, difference)==length(difference) #slope monotonically increases
            Sign = 2 #beneficial
        elseif count(i->i<=0, difference)==length(difference) #slope monotonically decreases
            Sign = -1 #deleterious
        else #if slope is not monotonic
            sflip = signflip1(difference) #for sign inversion, the derivative should only change sign once
            if sflip == 1
                Sign = 1 #sign inversion
            else
                Sign = 3 #unknown
            end
        end
    return(Sign)
end

function location(npfix)
    #classifies sign based on the location of points relative to the neutral benchmark
    if count(i->i==1, npfix)==length(npfix) #if all points == neutral
        Sign = 0 #neutral
    elseif count(i->i>=1, npfix)==length(npfix) #if all points >= neutral
        Sign = 2 #beneficial
    elseif count(i->i<=1, npfix)==length(npfix) #if all points <= neutral
        Sign = -1 #deleterious
    else #sign inversion means that all npfix less than some N will be deleterious, and all npfix greater than that N will be beneficial. therefore, the sign of selection should only "flip" twice
        sflip = signflip(npfix)
        if sflip == 2
            Sign = 1 #sign inversion
        else
            Sign = 3 #unknown
        end
    end
    return(Sign)
end

function neutral(npfix, allN, ci, reps)
    #determines if points are significantly different from 1
    npf = []
    for i in 1:length(npfix)
        bounds = ci[i]
        obs = npfix[i]/allN[i]*reps
        if bounds[1] < obs < bounds[2]
            push!(npf,1.0)
        else
            push!(npf,npfix[i])
        end
    end
    return(npf)
end

function interp1(npfix, allN)
    #interpolates the population size where the npfix curve crosses the neutral benchmark
    npfix = Interpolations.deduplicate_knots!(npfix)
    (~, idx) = findmin(npfix)
    interp = LinearInterpolation(npfix[idx:end],allN[idx:end])
    return (interp(1))
end

function signselect(npfix, allN, ci, reps)
    npf = neutral(npfix, allN, ci, reps) #first, evaluate if each point is significantly different from neutral
    m = monotonic(npf) #does the curve monotonically increase (beneficial) or monotonically decrease (deleterious)
    l = location(npf) #are all points >= 1 (beneficial) or <= 1 (deleterious)
    Ncrit = 0
    if m==3 || l==3 #if neither monotonic or location functions can classify, return "unknown sign"
        sign = 3
    elseif m==l #if both methods agree
        sign = m
        if sign == 1 #if the point exhibits sign inversion, interpolate ncrit
            Ncrit = interp1(npf,allN)
        else
            Ncrit = 0
        end
    else
        sign = 3 #if 2 methods disagree, return "unknown sign"
    end
    return(sign, Ncrit)
end

function confi(N, reps) #construct a neutral confidence interval at each population size to determine if an npfix value is significantly different from neutral
    b = Binomial(reps,1/N) #binomial distribution with p = probability of fixation for a neutral mutation
    p = pdf(b)
    c = []
    for i in 1:length(p)
       x = cdf(b,i)
       push!(c,x)
   end
   l = 1:length(c)
   cr = 1 .-c
   interp = LinearInterpolation(c,l)
   return([interp(0.01), interp(0.99)]) #0.01 confidence interval
end

day = string(now())

ci = [] #construct a neutral confidence interval at each population size to determine if an npfix value is significantly different from neutral
for N in allN
       c = confi(N,reps)
       push!(ci,c)
end

out = open("ncritspace_$etype.$ptype.$day.csv", "w") #create an output file
close(out)

out = open("ncritspace_npfix_$etype.$ptype.$day.csv", "w") #creates secondary output file that stores full NPFix vectors, in case the auto sorter cannot assign a sign of selection
close(out)

for pSpec in allpSpec
    for s in allS
        if pSpec == 1.0
            wc = s #if pSpec = 1, then wC will never be used, so a dummy value is assigned
        else
            wc = (5*pSpec - 2*(4*s^2 + (9*pSpec^2)/4)^(1/2))/(4*pSpec - 4) #calculates fitness of conservative phenotype necessary to yield GMF given pSpec
        end
        Npfix = Float64[]
        if wc > 0 #if there exists a possible wC for that point in paramspace. If not, the point is left blank
            for N in allN #repeats this process at each population size
                c = Atomic{Int64}(0) #counts number of replicates that reach fixation
                @threads for run = 1:reps
                    r = simulate(N, pA, pSpec, wc, envs, phenos) #c increases by 1 for each rep that reaches fixation
                    atomic_add!(c, Int(r))
                end
                push!(Npfix,Float64((c[]/reps)*N))
            end
            #after the full npfix curve is found, the sign of selection + ncrit are calculated
            (sign, Ncrit) = signselect(Npfix, allN, ci, reps)
        else #if there is no possible wC for this point in paramspace, the point is left blank
            sign = "NA"
            Ncrit = "NA"
        end

        out = open("ncritspace_$etype.$ptype.$day.csv", "a")
        write(out, join([pSpec, s, sign, Ncrit], ","), "\n")
        close(out)

        out = open("ncritspace_npfix_$etype.$ptype.$day.csv", "a")
        write(out, join([pSpec, wc, Npfix], ","), "\n")
        close(out)
    end
    println("$pSpec done")
end

tock()
