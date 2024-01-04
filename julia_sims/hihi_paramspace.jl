using Distributions
using Interpolations
using TickTock

tick()


filename = "hih_paramspace.csv"
println(filename)

pA = 0.5 #probability of being in env. A
allpSpec = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] #proportion of BH offspring with specialist phenotype
all_delta_gmf = [-0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075, 0.1. 0.125]

maxn = 6 #log10 maximum population sizes
a=range(0, stop=maxn, length = 20)
allN=[convert(Int64,floor(10^i)) for i in a] #creates vector of 20 population sizes from 1 to 10^6 evenly spaced on a log scale

out = open(filename, "w")
write(out, join(["model", "pSpec", "delta_gmf", "sign", "Ncrit", join(allN, ",")], ","), "\n") #creates output file with desired column names
close(out)

function beta_fitness(pSpec, gmf) #fitness of the b environment "beta" specialist in both environments
    #all_strat = ["hi/hi", "hi/lo"]
    if pSpec == 1 #if pSpec = 1, then no individuals will ever adopt the beta specialist phenotype
        if gmf == 1 #therefore, the only possible value of gmf = 0
            w_beta = [1, 1] #placeholder fitness values are assigned, but never used
        else
            w_beta = "skip" #the combination of pspec & gmf are impossible to achieve, so the point in param space is skipped
        end
    elseif pSpec == 0 #if pSPec == 0, then all individuals adopt the beta specialist, and it is no longer a bet hedger
        if gmf == 1
            w_beta = [0.5, 2] #pure beta specialist has opposite fitness of the alpha specialist
        else
            w_beta = "skip" #the combination of pspec & gmf are impossible to achieve, so the point in param space is skipped
        end
    else
        #first, we attempt to set w_beta_B = 2, the same high fitness val of the alpha specialist. we then back caulculate w_beta_B given all other params
        w_beta_A = (-(gmf^2)/(pSpec/2 - 2*(pSpec - 1)) + 2*pSpec)/(pSpec - 1)
        if w_beta_A >=0 #if w_beta_A is a non negative, and therefore possible, fitness value, we're done
            w_beta_B = 2
            w_beta = [w_beta_A, w_beta_B]
        else #but, if w_beta_A is negative, than w_beta_B = 2 is too high
            #we next try to set w_beta_A = 0.01, and back caulculate w_beta_B
            w_beta_B = (- (gmf^2)/(2*pSpec - 0.01*(pSpec - 1)) + pSpec/2)/(pSpec - 1)
            if w_beta_B >= 0 #if w_beta_B is a non negative, and therefore possible, fitness value, we're done
                w_beta_A = 0.01
                w_beta = [w_beta_A, w_beta_B]
            else #if w_beta_B is still negative, then we assume the combination of pspec & gmf are impossible to achieve, so the point in param space is skipped
                w_beta = "skip"
            end
        end
    end
    return w_beta
end

function env_select(pA) #randomly generates environment
    env = rand()
    return env <= pA
    #if returns true ==> in env. A; if returns false ==> in env. B
end

function pheno(pSpec, BHcount) #randomly generates phenotypic ddistribution of bet hedger offspring
    Alpha = rand(Binomial(BHcount, pSpec))
    #generates number of bet hedgers with specialist phenotype
    Beta = BHcount-Alpha
    return [Alpha, Beta]
    #returns a list with the # of bet hedgers with the alpha specialist and beta specialist phenotypes respectively
end

function average_fitness(PopNum, wBH, wWT)
    return (PopNum[1]*wBH+PopNum[2]*wWT)/(sum(PopNum))
end

function reproduction(env, PopNum, PopSize, spec, cons, phenotype_counts)
    #simulates random wright-fisher reproduction
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

function generation(PopNum, alpha, beta, pA, pSpec)
    env = env_select(pA) #determines the environment for this generation
    phenotype_counts = pheno(pSpec, PopNum[1]) #determines the realized phenotypic makeup of the BH
    PopNum = reproduction(env, PopNum, sum(PopNum), alpha, beta, phenotype_counts)
    #simulates random wright fisher reproduction
    return PopNum
end

function simulate(PopSize::Int64, pA::Float64, pSpec::Float64, w_beta_A::Float64, w_beta_B::Float64)
    InitNum = 1 #number of bh to start
    PopNum = [InitNum,PopSize-InitNum]
    #PopNum = [# of BH, # of WT]
    alpha = [2, 0.5] #alpha specialist phenotype in Env. A and Env. B respectively
    beta = [w_beta_A, w_beta_B] #beta specialist phenotype in Env. A and Env. B respectively
    generations::Int64 = 1
    while 0<PopNum[1]<PopSize
        #until the bet hedger reaches loss (count == 0) or fixation (count == PopSize)
        PopNum=generation(PopNum, alpha, beta, pA, pSpec)
        generations+=1
    end
    return PopNum[1]==PopSize
    #returns True if BH reaches fixation, False if BH is lost
end

function signflip1(vec) #for determining the sign of selection
    #calculates the number of times the [] changes sign
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

function signflip(vec) #for determining the sign of selection
    #calculates the number of times the slope changes sign
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

function monotonic(npfix) #determines if the npfix vector's slope is monotonic
    difference = [] #slope between each pair of points
    for i in collect(2:length(npfix))
        push!(difference, npfix[i]-npfix[i-1])
    end
        if count(i->i==0, difference)==length(difference) #if the slope is always == 0
            Sign = 0 #neutral
        elseif count(i->i>=0, difference)==length(difference) #if the slope is always positive, or increasing
            Sign = 2 #beneficial
        elseif count(i->i<=0, difference)==length(difference) #if the slope is always negative, or decreasing
            Sign = -1 #deleterious
        else
            sflip = signflip1(difference) #how many times does the slope change sign?
            if sflip == 1 #the slope should change sign only once if sign inversion occurs.
                Sign = 1 #sign inversion
            else
                Sign = 3 #unknown
            end
        end
    return(Sign)
end

function location(npfix)
    if count(i->i==1, npfix)==length(npfix) #if all points are not significantly different from neutral
        Sign = 0 #neutral
    elseif count(i->i>=1, npfix)==length(npfix) #if all points are greater than or equal to neutral
        Sign = 2 #beneficial
    elseif count(i->i<=1, npfix)==length(npfix) #if all points are less than or equal to neutral
        Sign = -1 #deleterious
    else
        sflip = signflip(npfix) #how many times does the location of points cross the neutral expectation?
        if sflip == 2 #should change twice --> once around N=1 (neutral --> del) and once around ncrit (del --> ben/neutral)
            Sign = 1 #sign inversion
        else
            Sign = 3 #unknown
        end
    end
    return(Sign)
end

function neutral(npfix, allN, ci) #all points are evaluated as being significantly different from the neutral expectation
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

function interp1(npfix, allN) #interpolate value of NCrit
    npfix = Interpolations.deduplicate_knots!(npfix)
    (~, idx) = findmin(npfix) #narrows down pop sizes
    interp = LinearInterpolation(npfix[idx:end],allN[idx:end])
    return (interp(1))
end

function signselect(npfix, allN, ci) #determine sign of selection for the npfix vector
    npf = neutral(npfix, allN, ci) #first, evaluate whether all points are different from neutral expectation
    m = monotonic(npf) #check sign of selection using the slope
    l = location(npf) #check the sign of selection using location of points
    Ncrit = 0 #assign empty value of Ncrit
    if m==3 || l==3 #if both monotonic and location functions assigned are "unknown"
        sign = 3
    elseif m==l #if both monotonic and location functions agree, return determined sign
        sign = m
        if sign == 1 #if the vector exhibits sign inversion, calculate Ncrit
            Ncrit = interp1(npf,allN)
        else
            Ncrit = 0
        end
    else #if the two functions disagree, assign as unknown
        sign = 3
    end
    return(sign, Ncrit)
end

function confi(N, reps) #creates a confidence interval for neutral expectation
    b = Binomial(reps,1/N) #neutral expectaion = 1/N
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

for pSpec in allpSpec
    for delta_gmf in all_delta_gmf
        gmf = 1+delta_gmf
        w_beta = beta_fitness(pSpec, gmf)
        if w_beta == "skip"
            println("Params skipped")
        else
            Npfix = Float64[]
            for N in allN #repeats this process at each population size
                if N < 10^4 #more reps == more time. this makes sure that you have just enough reps for significance
                    reps = 1000*10^4
                elseif N < 10^5
                    reps = 1000*10^5
                else
                    reps = 1000*10^6
                end
                c = 0
                for run = 1:reps
                    c += simulate(N, pA, pSpec, w_beta[1], w_beta[2])
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
    end
end

tock()
