using Distributions
using Interpolations
using TickTock

tick()

#day = string(now())

##Args --> $filename $SLURM_ARRAY_TASK_ID $pSpec $s
filename = join([ARGS[1], ARGS[2], "hihi", "csv"], ".")
println(filename)

pA = 0.5 #probability of being in env. A
#allpSpec = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] #proportion of BH offspring with specialist phenotype
pSpec = parse(Float64,ARGS[3])
#allS = [0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0]
s = parse(Float64,ARGS[4])

maxn = 6
a=range(0, stop=maxn, length = 15)
allN=[convert(Int64,floor(10^i)) for i in a]

out = open(filename, "w")
write(out, join(["model", "pSpec", "s", "sign", "Ncrit", join(allN, ",")], ","), "\n")
close(out)

function beta_fitness(pSpec, gmf, strat_model)
    #all_strat = ["hi/hi", "hi/lo"]
    if pSpec == 1
        if gmf == 1
            wc = [1, 1]
        else
            wc = "skip"
        end
    elseif strat_model == "hi/hi"
        wc_A = (-(gmf^2)/(pSpec/2 - 2*(pSpec - 1)) + 2*pSpec)/(pSpec - 1)
        if wc_A >=0
            wc_B = 2
            wc = [wc_A, wc_B]
        else
            wc_B = (- (gmf^2)/(2*pSpec - 0.01*(pSpec - 1)) + pSpec/2)/(pSpec - 1)
            if wc_B >= 0
                wc_A = 0.01
                wc = [wc_A, wc_B]
            else
                wc = "skip"
            end
        end
    else
        wcc = round((5*pSpec - 2*(4*gmf^2 + (9*pSpec^2)/4)^(1/2))/(4*pSpec - 4), digits = 5)
        if wc > 0
            wc = [wcc, wcc] #fitness in a, fitness in b
        else
            wc = "skip"
        end
    end
    return wc
end

function env_select(pA)
    env = rand()
    return env <= pA
    #if returns true ==> in env. A; if returns false ==> in env. B
end

function pheno(pSpec, BHcount)
    SpecPhe = rand(Binomial(BHcount, pSpec))
    #generates number of bet hedgers with specialist phenotype
    ConsPhe = BHcount-SpecPhe
    return [SpecPhe, ConsPhe]
    #returns a list with the # of bet hedgers with specialist and conservative phenotypes respectively
end

function average_fitness(PopNum, wBH, wWT)
    return (PopNum[1]*wBH+PopNum[2]*wWT)/(sum(PopNum))
end

function reproduction(env, PopNum, PopSize, spec, cons, phenotype_counts)
    wBH = Float64 #bet hedger fitness
    wWT = Float64 #wild type specialist fitness
    if env
        #if in env. A
        wBH = phenotype_counts[1]/PopNum[1]*spec[1]+phenotype_counts[2]/PopNum[1]*cons[1]
        wWT = spec[1]
    else
        #else in env. B
        wBH = phenotype_counts[1]/PopNum[1]*spec[2]+phenotype_counts[2]/PopNum[1]*cons[2]
        wWT = spec[2]
    end
    wBar = average_fitness(PopNum, wBH, wWT)
    p=PopNum[1]/PopSize * wBH/wBar
    PopNum[1] = rand(Binomial(PopSize, p))
    PopNum[2] = PopSize-PopNum[1]
    return PopNum
end

function generation(PopNum, spec, cons, pA, pSpec)
    env = env_select(pA) #determines the environment for this generation
    phenotype_counts = pheno(pSpec, PopNum[1]) #determines the realized phenotypic makeup of the BH
    PopNum = reproduction(env, PopNum, sum(PopNum), spec, cons, phenotype_counts)
    #simulates random wright fisher reproduction
    return PopNum
end

function simulate(PopSize::Int64, pA::Float64, pSpec::Float64, wc_A::Float64, wc_B::Float64)
    InitNum = 1 #number of bh to start
    PopNum = [InitNum,PopSize-InitNum]
    #PopNum = [# of BH, # of WT]
    spec = [2, 0.5] #specialist phenotype in Env. A and Env. B respectively
    cons = [wc_A, wc_B] #conservative phenotype in Env. A and Env. B respectively
    generations::Int64 = 1
    while 0<PopNum[1]<PopSize
        PopNum=generation(PopNum, spec, cons, pA, pSpec)
        generations+=1
    end
    return PopNum[1]==PopSize
    #returns True if BH win, False otherwise & number of generations elapsed
end

function signflip1(vec)
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

function signflip(vec)
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
    difference = []
    for i in collect(2:length(npfix))
        push!(difference, npfix[i]-npfix[i-1])
    end
        if count(i->i==0, difference)==length(difference)
            Sign = 0 #neutral
        elseif count(i->i>=0, difference)==length(difference)
            Sign = 2 #beneficial
        elseif count(i->i<=0, difference)==length(difference)
            Sign = -1 #deleterious
        else
            sflip = signflip1(difference)
            if sflip == 1
                Sign = 1 #sign inversion
            else
                Sign = 3 #unknown
            end
        end
    return(Sign)
end

function location(npfix)
    if count(i->i==1, npfix)==length(npfix)
        Sign = 0 #neutral
    elseif count(i->i>=1, npfix)==length(npfix)
        Sign = 2 #beneficial
    elseif count(i->i<=1, npfix)==length(npfix)
        Sign = -1 #deleterious
    else
        sflip = signflip(npfix)
        if sflip == 2
            Sign = 1 #sign inversion
        else
            Sign = 3 #unknown
        end
    end
    return(Sign)
end

function neutral(npfix, allN, ci)
    npf = []
    for i in 1:length(npfix)
        bounds = ci[i]
        if bounds[1] < npfix[i] < bounds[2]
            push!(npf,1.0)
        else
            push!(npf,npfix[i])
        end
    end
    return(npf)
end

function interp1(npfix, allN)
    npfix = Interpolations.deduplicate_knots!(npfix)
    (~, idx) = findmin(npfix)
    interp = LinearInterpolation(npfix[idx:end],allN[idx:end])
    return (interp(1))
end

function signselect(npfix, allN, ci)
    npf = neutral(npfix, allN, ci)
    m = monotonic(npf)
    l = location(npf)
    Ncrit = 0
    if m==3 || l==3
        sign = 3
    elseif m==l
        sign = m
        if sign == 1
            Ncrit = interp1(npf,allN)
        else
            Ncrit = 0
        end
    else
        sign = 3
    end
    return(sign, Ncrit)
end

function confi(N, reps)
    b = Binomial(reps,1/N)
    p = pdf(b)
    c = []
    for i in 1:length(p)
       x = cdf(b,i)
       push!(c,x)
   end
   l = 1:length(c)
   cr = 1 .-c
   interp = LinearInterpolation(c,l)
   return([(interp(0.01)/reps)*N, (interp(0.99)/reps)*N])
end

ci = []

gmf = 1+s
wC = beta_fitness(pSpec, gmf, "hi/hi")
if wC == "skip"
    println("Params skipped")
else
    Npfix = Float64[]
    for N in allN #repeats this process at each population size
        reps = 10000
        if N < 10^4
            reps = 1000*10^4
        elseif N < 10^5
            reps = 1000*10^5
        else
            reps = 1000*10^6
        end
        c = 0
        for run = 1:reps
            c += simulate(N, pA, pSpec, wC[1], wC[2])
        end
        push!(Npfix,Float64((c/reps)*N))
        conf_int = confi(N,reps)
        push!(ci,conf_int)
    end

    (sign_sel, Ncrit) = signselect(Npfix, allN, ci)

    out = open(filename, "a")
    write(out, join(["hi/hi", pSpec, s, sign_sel, Ncrit, join(Npfix, ",")], ","), "\n")
    close(out)
end

tock()