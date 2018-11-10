var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Cuba",
    "title": "Cuba",
    "category": "page",
    "text": ""
},

{
    "location": "#Cuba-1",
    "page": "Cuba",
    "title": "Cuba",
    "category": "section",
    "text": "DocTestSetup = quote\n    using Cuba\nend"
},

{
    "location": "#Introduction-1",
    "page": "Cuba",
    "title": "Introduction",
    "category": "section",
    "text": "Cuba.jl is a Julia library for multidimensional numerical integration of real-valued functions of real arguments, using different algorithms.This is just a Julia wrapper around the C Cuba library, version 4.2, by Thomas Hahn.  All the credits goes to him for the underlying functions, blame me for any problem with the Julia interface.All algorithms provided by Cuba library are supported in Cuba.jl:Vegas:\nBasic integration method Type Variance reduction\nSobol quasi-random sample Monte Carlo importance sampling\nMersenne Twister pseudo-random sample \" \"\nRanlux pseudo-random sample \" \"\nSuave\nBasic integration method Type Variance reduction\nSobol quasi-random sample Monte Carlo globally adaptive subdivision and importance sampling\nMersenne Twister pseudo-random sample \" \"\nRanlux pseudo-random sample \" \"\nDivonne\nBasic integration method Type Variance reduction\nKorobov quasi-random sample Monte Carlo stratified sampling aided by methods from numerical optimization\nSobol quasi-random sample \" \"\nMersenne Twister pseudo-random sample \" \"\nRanlux pseudo-random sample \" \"\ncubature rules deterministic \"\nCuhre\nBasic integration method Type Variance reduction\ncubature rules deterministic globally adaptive subdivisionFor more details on the algorithms see the manual included in Cuba library and available in deps/usr/share/cuba.pdf after successful installation of Cuba.jl.Integration is always performed on the n-dimensional unit hypercube 0 1^n.tip: Tip\nIf you want to compute an integral over a different set, you have to scale the integrand function in order to have an equivalent integral on 0 1^n using substitution rules. For example, recall that in one dimensionint_a^b f(x)mathrmdx = int_0^1 f(a + (b - a) y) (b -\na)mathrmdywhere the final (b - a) is the one-dimensional version of the Jacobian.Integration over a semi-infinite or an inifite domain is a bit trickier, but you can follow this advice from Steven G. Johnson: to compute an integral over a semi-infinite interval, you can perform the change of variables x=a+y(1-y):int_a^infty f(x)mathrmdx = int_0^1\nfleft(a + fracy1 - yright)frac1(1 - y)^2mathrmdyFor an infinite interval, you can perform the change of variables x=(2y - 1)((1 - y)y):int_-infty^infty f(x)mathrmdx = int_0^1\nfleft(frac2y - 1(1 - y)yright)frac2y^2 - 2y + 1(1 -\ny)^2y^2mathrmdyIn addition, recall that for an even function int_-infty^infty f(x)mathrmdx = 2int_0^inftyf(x)mathrmdx, while the integral of an odd function over the infinite interval (-infty infty) is zero.All this generalizes straightforwardly to more than one dimension. In Examples section you can find the computation of a 3-dimensional integral with non-constant boundaries using Cuba.jl and two integrals over infinite domains.Cuba.jl is available for GNU/Linux, FreeBSD, Mac OS, and Windows (i686 and x86_64 architectures)."
},

{
    "location": "#Installation-1",
    "page": "Cuba",
    "title": "Installation",
    "category": "section",
    "text": "Cuba.jl is available for Julia 0.7 and later versions, and can be installed with Julia built-in package manager. In a Julia session run the commandsjulia> Pkg.update()\njulia> Pkg.add(\"Cuba\")Older versions are also available for Julia 0.4-0.6."
},

{
    "location": "#Cuba.vegas",
    "page": "Cuba",
    "title": "Cuba.vegas",
    "category": "function",
    "text": "vegas(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions\n\nCalculate integral of integrand over the unit hypercube in ndim dimensions using Vegas algorithm.  integrand is a vectorial function with ncomp components.  ndim and ncomp default to 1.\n\nAccepted keywords:\n\nnvec\nrtol\natol\nflags\nseed\nminevals\nmaxevals\nnstart\nnincrease\nnbatch\ngridno\nstatefile\nspin\n\n\n\n\n\n"
},

{
    "location": "#Cuba.suave",
    "page": "Cuba",
    "title": "Cuba.suave",
    "category": "function",
    "text": "suave(integrand, ndim=1, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions\n\nCalculate integral of integrand over the unit hypercube in ndim dimensions using Suave algorithm.  integrand is a vectorial function with ncomp components. ndim and ncomp default to 1.\n\nAccepted keywords:\n\nnvec\nrtol\natol\nflags\nseed\nminevals\nmaxevals\nnnew\nnmin\nflatness\nstatefile\nspin\n\n\n\n\n\n"
},

{
    "location": "#Cuba.divonne",
    "page": "Cuba",
    "title": "Cuba.divonne",
    "category": "function",
    "text": "divonne(integrand, ndim=2, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions\n\nCalculate integral of integrand over the unit hypercube in ndim dimensions using Divonne algorithm.  integrand is a vectorial function with ncomp components. ncomp defaults to 1, ndim defaults to 2 and must be ≥ 2.\n\nAccepted keywords:\n\nnvec\nrtol\natol\nflags\nseed\nminevals\nmaxevals\nkey1\nkey2\nkey3\nmaxpass\nborder\nmaxchisq\nmindeviation\nngiven\nldxgiven\nxgiven\nnextra\npeakfinder\nstatefile\nspin\n\n\n\n\n\n"
},

{
    "location": "#Cuba.cuhre",
    "page": "Cuba",
    "title": "Cuba.cuhre",
    "category": "function",
    "text": "cuhre(integrand, ndim=2, ncomp=1[, keywords]) -> integral, error, probability, neval, fail, nregions\n\nCalculate integral of integrand over the unit hypercube in ndim dimensions using Cuhre algorithm.  integrand is a vectorial function with ncomp components.  ncomp defaults to 1, ndim defaults to 2 and must be ≥ 2.\n\nAccepted keywords:\n\nnvec\nrtol\natol\nflags\nminevals\nmaxevals\nkey\nstatefile\nspin\n\n\n\n\n\n"
},

{
    "location": "#Usage-1",
    "page": "Cuba",
    "title": "Usage",
    "category": "section",
    "text": "After installing the package, runusing Cubaor put this command into your Julia script.Cuba.jl provides the following functions to integrate:vegas\nsuave\ndivonne\ncuhreLarge parts of the following sections are borrowed from Cuba manual. Refer to it for more information on the details.Cuba.jl wraps the 64-bit integers functions of Cuba library, in order to push the range of certain counters to its full extent. In detail, the following arguments:for Vegas: nvec, minevals, maxevals, nstart, nincrease,   nbatch, neval,\nfor Suave: nvec, minevals, maxevals, nnew, nmin, neval,\nfor Divonne: nvec, minevals, maxevals, ngiven, nextra,   neval,\nfor Cuhre: nvec, minevals, maxevals, neval,are passed to the Cuba library as 64-bit integers, so they are limited to be at mostjulia> typemax(Int64)\n9223372036854775807There is no way to overcome this limit. See the following sections for the meaning of each argument."
},

{
    "location": "#Arguments-1",
    "page": "Cuba",
    "title": "Arguments",
    "category": "section",
    "text": "The only mandatory argument of integrator functions is:integrand (type: Function): the function to be integratedOptional positional arguments are:ndim (type: Integer): the number of dimensions of the   integratation domain. If omitted, defaults to 1 in vegas and   suave, to 2 in divonne and cuhre. Note: ndim must be at   least 2 with the latest two methods.\nncomp (type: Integer): the number of components of the   integrand. Default to 1 if omittedintegrand should be a function integrand(x, f) taking two arguments:the input vector x of length ndim\nthe output vector f of length ncomp, used to set the value of   each component of the integrand at point xx and f are matrices with dimensions (ndim, nvec) and (ncomp, nvec), respectively, when nvec > 1. See the Vectorization section below for more information.Also anonymous functions are allowed as integrand. For those familiar with Cubature.jl package, this is the same syntax used for integrating vector-valued functions.For example, the integralint_0^1 cos (x) mathrmdx = sin(1) = 08414709848078965can be computed with one of the following commandsjulia> vegas((x, f) -> f[1] = cos(x[1]))\nComponent:\n 1: 0.8414910005259609 ± 5.2708169787733e-5 (prob.: 0.028607201257039333)\nIntegrand evaluations: 13500\nNumber of subregions:  0\nNote: The desired accuracy was reached\n\njulia> suave((x, f) -> f[1] = cos(x[1]))\nComponent:\n 1: 0.8411523690658836 ± 8.357995611133613e-5 (prob.: 1.0)\nIntegrand evaluations: 22000\nNumber of subregions:  22\nNote: The desired accuracy was reached\n\njulia> divonne((x, f) -> f[1] = cos(x[1]))\nComponent:\n 1: 0.841468071955942 ± 5.3955070531551656e-5 (prob.: 0.0)\nIntegrand evaluations: 1686\nNumber of subregions:  14\nNote: The desired accuracy was reached\n\njulia> cuhre((x, f) -> f[1] = cos(x[1]))\nComponent:\n 1: 0.8414709848078966 ± 2.2204460420128823e-16 (prob.: 3.443539937576958e-5)\nIntegrand evaluations: 195\nNumber of subregions:  2\nNote: The desired accuracy was reachedIn section Examples you can find more complete examples.  Note that x and f are both arrays with type Float64, so Cuba.jl can be used to integrate real-valued functions of real arguments. See how to work with a complex integrand.Note: if you used Cuba.jl until version 0.0.4, be aware that the user interface has been reworked in version 0.0.5 in a backward incompatible way."
},

{
    "location": "#Optional-Keywords-1",
    "page": "Cuba",
    "title": "Optional Keywords",
    "category": "section",
    "text": "All other arguments required by Cuba integrator routines can be passed as optional keywords. Cuba.jl uses some reasonable default values in order to enable users to invoke integrator functions with a minimal set of arguments. Anyway, if you want to make sure future changes to some default values of keywords will not affect your current script, explicitely specify the value of the keywords."
},

{
    "location": "#Common-Keywords-1",
    "page": "Cuba",
    "title": "Common Keywords",
    "category": "section",
    "text": "These are optional keywords common to all functions:nvec (type: Integer, default: 1): the maximum number of points to be   given to the integrand routine in each invocation. Usually this is 1 but if   the integrand can profit from e.g. Single Instruction Multiple Data (SIMD)   vectorization, a larger value can be chosen. See Vectorization   section.\nrtol (type: Real, default: 1e-4), and atol (type: Real,   default: 1e-12): the requested relative   (varepsilon_textrel) and absolute   (varepsilon_textabs) accuracies. The integrator tries to   find an estimate hatI for the integral I which for every   component c fulfills hatI_c - I_cleq   max(varepsilon_textabs varepsilon_textrel I_c).\nflags (type: Integer, default: 0): flags governing the   integration:\nBits 0 and 1 are taken as the verbosity level, i.e. 0 to 3,   unless the CUBAVERBOSE environment variable contains an even   higher value (used for debugging).\nLevel 0 does not print any output, level 1 prints \"reasonable\"   information on the progress of the integration, level 2 also echoes   the input parameters, and level 3 further prints the subregion results   (if applicable).\nBit 2 = 0: all sets of samples collected on a subregion during   the various iterations or phases contribute to the final result.\nBit 2 = 1, only the last (largest) set of samples is used in   the final result.\n(Vegas and Suave only)\nBit 3 = 0, apply additional smoothing to the importance   function, this moderately improves convergence for many   integrands.\nBit 3 = 1, use the importance function without smoothing, this   should be chosen if the integrand has sharp edges.\nBit 4 = 0, delete the state file (if one is chosen) when the   integration terminates successfully.\nBit 4 = 1, retain the state file.\n(Vegas only)\nBit 5 = 0, take the integrator\'s state from the state file,   if one is present.\nBit 5 = 1, reset the integrator\'s state even if a state file   is present, i.e. keep only the grid. Together with Bit 4 this   allows a grid adapted by one integration to be used for another   integrand.\nBits 8–31 =: level determines the random-number generator.\nTo select e.g. last samples only and verbosity level 2, pass   6 = 4 + 2 for the flags.\nseed (type: Integer, default: 0): the seed for the pseudo-random-number generator. This keyword is not available for cuhre. The random-number generator is chosen as follows:\nseed level (bits 8–31 of flags) Generator\nzero N/A Sobol (quasi-random)\nnon-zero zero Mersenne Twister (pseudo-random)\nnon-zero non-zero Ranlux (pseudo-random)\nRanlux implements Marsaglia and Zaman\'s 24-bit RCARRY algorithm with generation period p, i.e. for every 24 generated numbers used, another p - 24 are skipped. The luxury level is encoded in level as follows:\nLevel 1 (p = 48): very long period, passes the gap test but fails spectral test.\nLevel 2 (p = 97): passes all known tests, but theoretically still defective.\nLevel 3 (p = 223): any theoretically possible correlations have very small chance of being observed.\nLevel 4 (p = 389): highest possible luxury, all 24 bits chaotic.\nLevels 5–23 default to 3, values above 24 directly specify the period p. Note that Ranlux\'s original level 0, (mis)used for selecting Mersenne Twister in Cuba, is equivalent to level = 24.\nminevals (type: Real, default: 0): the minimum number of   integrand evaluations required\nmaxevals (type: Real, default: 1000000): the (approximate)   maximum number of integrand evaluations allowed\nstatefile (type: AbstractString, default: \"\"): a filename for   storing the internal state. To not store the internal state, put   \"\" (empty string, this is the default) or C_NULL (C null   pointer).\nCuba can store its entire internal state (i.e. all the information   to resume an interrupted integration) in an external file. The state   file is updated after every iteration. If, on a subsequent   invocation, a Cuba routine finds a file of the specified name, it   loads the internal state and continues from the point it left off.   Needless to say, using an existing state file with a different   integrand generally leads to wrong results.\nThis feature is useful mainly to define \"check-points\" in   long-running integrations from which the calculation can be   restarted.\nOnce the integration reaches the prescribed accuracy, the state file   is removed, unless bit 4 of flags (see above) explicitly requests   that it be kept.\nspin (type: Ptr{Void}, default: C_NULL): this is the   placeholder for the \"spinning cores\" pointer. Cuba.jl does not   support parallelization, so this keyword should not have a value   different from C_NULL."
},

{
    "location": "#Vegas-Specific-Keywords-1",
    "page": "Cuba",
    "title": "Vegas-Specific Keywords",
    "category": "section",
    "text": "These optional keywords can be passed only to vegas:nstart (type: Integer, default: 1000): the number of integrand   evaluations per iteration to start with\nnincrease (type: Integer, default: 500): the increase in the   number of integrand evaluations per iteration\nnbatch (type: Integer, default: 1000): the batch size for   sampling\nVegas samples points not all at once, but in batches of size   nbatch, to avoid excessive memory consumption. 1000 is a   reasonable value, though it should not affect performance too much\ngridno (type: Integer, default: 0): the slot in the internal   grid table.\nIt may accelerate convergence to keep the grid accumulated during   one integration for the next one, if the integrands are reasonably   similar to each other. Vegas maintains an internal table with space   for ten grids for this purpose. The slot in this grid is specified   by gridno.\nIf a grid number between 1 and 10 is selected, the grid is not   discarded at the end of the integration, but stored in the   respective slot of the table for a future invocation. The grid is   only re-used if the dimension of the subsequent integration is the   same as the one it originates from.\nIn repeated invocations it may become necessary to flush a slot in   memory, in which case the negative of the grid number should be set."
},

{
    "location": "#Suave-Specific-Keywords-1",
    "page": "Cuba",
    "title": "Suave-Specific Keywords",
    "category": "section",
    "text": "These optional keywords can be passed only to suave:nnew (type: Integer, default: 1000): the number of new   integrand evaluations in each subdivision\nnmin (type: Integer, default: 2): the minimum number of   samples a former pass must contribute to a subregion to be   considered in that region\'s compound integral value. Increasing   nmin may reduce jumps in the chi^2 value\nflatness (type: Real, default: .25): the type of norm used to   compute the fluctuation of a sample. This determines how prominently   \"outliers\", i.e. individual samples with a large fluctuation,   figure in the total fluctuation, which in turn determines how a   region is split up. As suggested by its name, flatness should be   chosen large for \"flat\" integrands and small for \"volatile\"   integrands with high peaks. Note that since flatness appears in   the exponent, one should not use too large values (say, no more than   a few hundred) lest terms be truncated internally to prevent   overflow."
},

{
    "location": "#Divonne-Specific-Keywords-1",
    "page": "Cuba",
    "title": "Divonne-Specific Keywords",
    "category": "section",
    "text": "These optional keywords can be passed only to divonne:key1 (type: Integer, default: 47): determines sampling in the   partitioning phase: key1 = 7 9 11 13 selects the cubature   rule of degree key1. Note that the degree-11 rule is available   only in 3 dimensions, the degree-13 rule only in 2 dimensions.\nFor other values of key1, a quasi-random sample of n_1 =   verbkey1 points is used, where the sign of key1 determines   the type of sample,\nkey1  0, use a Korobov quasi-random sample,\nkey1  0, use a \"standard\" sample (a Sobol quasi-random   sample if seed = 0, otherwise a pseudo-random sample).\nkey2 (type: Integer, default: 1): determines sampling in   the final integration phase:\nkey2 = 7 9 11 13 selects the cubature rule of degree   key2. Note that the degree-11 rule is available only in 3   dimensions, the degree-13 rule only in 2 dimensions.\nFor other values of key2, a quasi-random sample is used, where   the sign of key2 determines the type of sample,\nkey2  0, use a Korobov quasi-random sample,\nkey2  0, use a \"standard\" sample (see description of   key1 above),\nand n_2 = verbkey2 determines the number of points,\nn_2geq 40, sample n_2 points,\nn_2  40, sample n_2n_textneed points, where   n_textneed is the number of points needed to reach   the prescribed accuracy, as estimated by Divonne from the   results of the partitioning phase\nkey3 (type: Integer, default: 1): sets the strategy for the   refinement phase:\nkey3 = 0, do not treat the subregion any further.\nkey3 = 1, split the subregion up once more.\nOtherwise, the subregion is sampled a third time with key3   specifying the sampling parameters exactly as key2 above.\nmaxpass (type: Integer, default: 5): controls the thoroughness   of the partitioning phase: The partitioning phase terminates when   the estimated total number of integrand evaluations (partitioning   plus final integration) does not decrease for maxpass successive   iterations.\nA decrease in points generally indicates that Divonne discovered new   structures of the integrand and was able to find a more effective   partitioning. maxpass can be understood as the number of   \"safety\" iterations that are performed before the partition is   accepted as final and counting consequently restarts at zero   whenever new structures are found.\nborder (type: Real, default: 0.): the width of the border of   the integration region. Points falling into this border region will   not be sampled directly, but will be extrapolated from two samples   from the interior. Use a non-zero border if the integrand function   cannot produce values directly on the integration boundary\nmaxchisq (type: Real, default: 10.): the chi^2 value a   single subregion is allowed to have in the final integration phase.   Regions which fail this chi^2 test and whose sample averages   differ by more than mindeviation move on to the refinement phase.\nmindeviation (type: Real, default: 0.25): a bound, given as   the fraction of the requested error of the entire integral, which   determines whether it is worthwhile further examining a region that   failed the chi^2 test. Only if the two sampling averages obtained   for the region differ by more than this bound is the region further   treated.\nngiven (type: Integer, default: 0): the number of points in   the xgiven array\nldxgiven (type: Integer, default: 0): the leading dimension of   xgiven, i.e. the offset between one point and the next in memory\nxgiven (type: AbstractArray{Real}, default:   zeros(Cdouble, ldxgiven, ngiven)): a list of points where the   integrand might have peaks. Divonne will consider these points when   partitioning the integration region. The idea here is to help the   integrator find the extrema of the integrand in the presence of very   narrow peaks. Even if only the approximate location of such peaks is   known, this can considerably speed up convergence.\nnextra (type: Integer, default: 0): the maximum number of   extra points the peak-finder subroutine will return. If nextra is   zero, peakfinder is not called and an arbitrary object may be   passed in its place, e.g. just 0\npeakfinder (type: Ptr{Void}, default: C_NULL): the peak-finder   subroutine"
},

{
    "location": "#Cuhre-Specific-Keyword-1",
    "page": "Cuba",
    "title": "Cuhre-Specific Keyword",
    "category": "section",
    "text": "This optional keyword can be passed only to cuhre:key (type: Integer, default: 0): chooses the basic integration   rule:\nkey = 7 9 11 13 selects the cubature rule of degree key.   Note that the degree-11 rule is available only in 3 dimensions, the   degree-13 rule only in 2 dimensions.\nFor other values, the default rule is taken, which is the degree-13   rule in 2 dimensions, the degree-11 rule in 3 dimensions, and the   degree-9 rule otherwise."
},

{
    "location": "#Output-1",
    "page": "Cuba",
    "title": "Output",
    "category": "section",
    "text": "The integrating functions vegas, suave, divonne, and cuhre return an Integral object whose fields areintegral :: Vector{Float64}\nerror    :: Vector{Float64}\nprobl    :: Vector{Float64}\nneval    :: Int64\nfail     :: Int32\nnregions :: Int32The first three fields are arrays with length ncomp, the last three ones are scalars. The Integral object can also be iterated over like a tuple. In particular, if you assign the output of integrator functions to the variable named result, you can access the value of the i-th component of the integral with result[1][i] or result.integral[i] and the associated error with result[2][i] or result.error[i].integral (type: Vector{Float64}, with ncomp components): the   integral of integrand over the unit hypercube\nerror (type: Vector{Float64}, with ncomp components): the   presumed absolute error for each component of integral\nprobability (type: Vector{Float64}, with ncomp components):   the chi^2 -probability (not the chi^2 -value itself!) that   error is not a reliable estimate of the true integration error. To   judge the reliability of the result expressed through prob,   remember that it is the null hypothesis that is tested by the   chi^2 test, which is that error is a reliable estimate. In   statistics, the null hypothesis may be rejected only if prob is   fairly close to unity, say prob 95\nneval (type: Int64): the actual number of integrand evaluations   needed\nfail (type: Int32): an error flag:\nfail = 0, the desired accuracy was reached\nfail = -1, dimension out of range\nfail > 0, the accuracy goal was not met within the allowed   maximum number of integrand evaluations. While Vegas, Suave, and   Cuhre simply return 1, Divonne can estimate the number of   points by which maxevals needs to be increased to reach the   desired accuracy and returns this value.\nnregions (type: Int32): the actual number of subregions needed   (always 0 in vegas)"
},

{
    "location": "#Vectorization-1",
    "page": "Cuba",
    "title": "Vectorization",
    "category": "section",
    "text": "Vectorization means evaluating the integrand function for several points at once. This is also known as Single Instruction Multiple Data (SIMD) paradigm and is different from ordinary parallelization where independent threads are executed concurrently. It is usually possible to employ vectorization on top of parallelization.Cuba.jl cannot automatically vectorize the integrand function, of course, but it does pass (up to) nvec points per integrand call (Common Keywords). This value need not correspond to the hardware vector length –computing several points in one call can also make sense e.g. if the computations have significant intermediate results in common.When nvec > 1, the input x is a matrix of dimensions (ndim, nvec), while the output f is a matrix with dimensions (ncomp, nvec). Vectorization can be used to evaluate more quickly the integrand function, for example by exploiting parallelism, thus speeding up computation of the integral. See the section Vectorized Function below for an example of a vectorized funcion.note: Disambiguation\nThe nbatch argument of vegas is related in purpose but not identical to nvec. It internally partitions the sampling done by Vegas but has no bearing on the number of points given to the integrand. On the other hand, it it pointless to choose nvec > nbatch for Vegas."
},

{
    "location": "#Examples-1",
    "page": "Cuba",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "#One-dimensional-integral-1",
    "page": "Cuba",
    "title": "One dimensional integral",
    "category": "section",
    "text": "The integrand ofint_0^1 fraclog(x)sqrtx mathrmdxhas an algebraic-logarithmic divergence for x = 0, but the integral is convergent and its value is -4. Cuba.jl integrator routines can handle this class of functions and you can easily compute the numerical approximation of this integral using one of the following commands:julia> vegas( (x,f) -> f[1] = log(x[1])/sqrt(x[1]))\nComponent:\n 1: -3.9981623937128465 ± 0.00044066437168409464 (prob.: 0.2843052968819913)\nIntegrand evaluations: 1007500\nNumber of subregions:  0\nNote: The accuracy was not met within the maximum number of evaluations\n\njulia> suave( (x,f) -> f[1] = log(x[1])/sqrt(x[1]))\nComponent:\n 1: -3.999976286717141 ± 0.00039504866661339003 (prob.: 1.0)\nIntegrand evaluations: 51000\nNumber of subregions:  51\nNote: The desired accuracy was reached\n\njulia> divonne( (x,f) -> f[1] = log(x[1])/sqrt(x[1]), atol = 1e-8, rtol = 1e-8)\nComponent:\n 1: -3.999999899620808 ± 2.1865962888459237e-7 (prob.: 0.0)\nIntegrand evaluations: 1002059\nNumber of subregions:  1582\nNote: The accuracy was not met within the maximum number of evaluations\nHint: Try increasing `maxevals` to 4884287\n\njulia> cuhre( (x,f) -> f[1] = log(x[1])/sqrt(x[1]))\nComponent:\n 1: -4.000000355067187 ± 0.0003395484028626406 (prob.: 0.0)\nIntegrand evaluations: 5915\nNumber of subregions:  46\nNote: The desired accuracy was reached"
},

{
    "location": "#Vector-valued-integrand-1",
    "page": "Cuba",
    "title": "Vector-valued integrand",
    "category": "section",
    "text": "Consider the integralintlimits_Omega\nboldsymbolf(xyz)mathrmdxmathrmdymathrmdzwhere Omega = 0 1^3 andboldsymbolf(xyz) = left(sin(x)cos(y)exp(z) exp(-(x^2 + y^2 +\nz^2)) frac11 - xyzright)In this case it is more convenient to write a simple Julia script to compute the above integraljulia> using Cuba, SpecialFunctions\n\njulia> function integrand(x, f)\n           f[1] = sin(x[1])*cos(x[2])*exp(x[3])\n           f[2] = exp(-(x[1]^2 + x[2]^2 + x[3]^2))\n           f[3] = 1/(1 - prod(x))\n       end\nintegrand (generic function with 1 method)\n\njulia> result, err = cuhre(integrand, 3, 3, atol=1e-12, rtol=1e-10);\n\njulia> answer = ((ℯ-1)*(1-cos(1))*sin(1), (sqrt(pi)*erf(1)/2)^3, zeta(3));\n\njulia> for i = 1:3\n           println(\"Component \", i)\n           println(\" Result of Cuba: \", result[i], \" ± \", err[i])\n           println(\" Exact result:   \", answer[i])\n           println(\" Actual error:   \", abs(result[i] - answer[i]))\n       end\nComponent 1\n Result of Cuba: 0.6646696797813745 ± 1.0056262721114345e-13\n Exact result:   0.6646696797813771\n Actual error:   2.6645352591003757e-15\nComponent 2\n Result of Cuba: 0.41653838588064585 ± 2.932867102879894e-11\n Exact result:   0.41653838588663805\n Actual error:   5.992206730809357e-12\nComponent 3\n Result of Cuba: 1.2020569031649704 ± 1.1958521782293645e-10\n Exact result:   1.2020569031595951\n Actual error:   5.375255796025158e-12"
},

{
    "location": "#Integral-with-non-constant-boundaries-1",
    "page": "Cuba",
    "title": "Integral with non-constant boundaries",
    "category": "section",
    "text": "The integralint_-y^yint_0^zint_0^pi\ncos(x)sin(y)exp(z)mathrmdxmathrmdymathrmdzhas non-constant boundaries. By applying the substitution rule repeatedly, you can scale the integrand function and get this equivalent integral over the fixed domain Omega = 0 1^3intlimits_Omega 2pi^3yz^2 cos(pi yz(2x - 1)) sin(pi yz)\nexp(pi z)mathrmdxmathrmdymathrmdzthat can be computed with Cuba.jl using the following Julia scriptjulia> using Cuba\n\njulia> function integrand(x, f)\n           f[1] = 2pi^3*x[2]*x[3]^2*cos(pi*x[2]*x[3]*(2*x[1] - 1.0))*\n                  sin(pi*x[2]*x[3])*exp(pi*x[3])\n       end\nintegrand (generic function with 1 method)\n\njulia> result, err = cuhre(integrand, 3, 1, atol=1e-12, rtol=1e-10);\n\njulia> answer = pi*ℯ^pi - (4ℯ^pi - 4)/5;\n\njulia> begin\n               println(\"Result of Cuba: \", result[1], \" ± \", err[1])\n               println(\"Exact result:   \", answer)\n               println(\"Actual error:   \", abs(result[1] - answer))\n       end\nResult of Cuba: 54.98607586826155 ± 5.4606062189698135e-9\nExact result:   54.98607586789537\nActual error:   3.6617819887396763e-10"
},

{
    "location": "#Integrals-over-Infinite-Domains-1",
    "page": "Cuba",
    "title": "Integrals over Infinite Domains",
    "category": "section",
    "text": "Cuba.jl assumes always as integration domain the hypercube 0 1^n, but we have seen that using integration by substitution we can calculate integrals over different domains as well. In the Introduction we also proposed two useful substitutions that can be employed to change an infinite or semi-infinite domain into a finite one.As a first example, consider the following integral with a semi-infinite domain:int_0^inftyfraclog(1 + x^2)1 + x^2mathrmdxwhose exact result is pilog 2. This can be computed as follows:julia> using Cuba\n\njulia> # The function we want to integrate over [0, ∞).\n\njulia> func(x) = log(1 + x^2)/(1 + x^2)\nfunc (generic function with 1 method)\n\njulia> # Scale the function in order to integrate over [0, 1].\n\njulia> function integrand(x, f)\n           f[1] = func(x[1]/(1 - x[1]))/(1 - x[1])^2\n       end\nintegrand (generic function with 1 method)\n\njulia> result, err = cuhre(integrand, atol = 1e-12, rtol = 1e-10);\n\njulia> answer = pi*log(2);\n\njulia> begin\n               println(\"Result of Cuba: \", result[1], \" ± \", err[1])\n               println(\"Exact result:   \", answer)\n               println(\"Actual error:   \", abs(result[1] - answer))\n       end\nResult of Cuba: 2.1775860903056885 ± 2.150398850102772e-10\nExact result:   2.177586090303602\nActual error:   2.086331107875594e-12Now we want to calculate this integral, over an infinite domainint_-infty^infty frac1 - cos xx^2mathrmdxwhich gives pi. You can calculate the result with the code below. Note that integrand function has value 12 for x=0, but you have to inform Julia about this.julia> using Cuba\n\njulia> # The function we want to integrate over (-∞, ∞).\n\njulia> func(x) = x==0 ? 0.5*one(x) : (1 - cos(x))/x^2\nfunc (generic function with 1 method)\n\njulia> # Scale the function in order to integrate over [0, 1].\n\njulia> function integrand(x, f)\n           f[1] = func((2*x[1] - 1)/x[1]/(1 - x[1])) *\n                   (2*x[1]^2 - 2*x[1] + 1)/x[1]^2/(1 - x[1])^2\n       end\nintegrand (generic function with 1 method)\n\njulia> result, err = cuhre(integrand, atol = 1e-7, rtol = 1e-7);\n\njulia> answer = float(pi);\n\njulia> begin\n               println(\"Result of Cuba: \", result[1], \" ± \", err[1])\n               println(\"Exact result:   \", answer)\n               println(\"Actual error:   \", abs(result[1] - answer))\n       end\nResult of Cuba: 3.1415928900555046 ± 2.050669142074607e-6\nExact result:   3.141592653589793\nActual error:   2.3646571145619077e-7"
},

{
    "location": "#Complex-integrand-1",
    "page": "Cuba",
    "title": "Complex integrand",
    "category": "section",
    "text": "As already explained, Cuba.jl operates on real quantities, so if you want to integrate a complex-valued function of complex arguments you have to treat complex quantities as 2-component arrays of real numbers.  For example, if you do not remember Euler\'s formula, you can compute this simple integralint_0^pi2 exp(mathrmi x)mathrmdxwith the following codejulia> using Cuba\n\njulia> function integrand(x, f)\n           # Complex integrand, scaled to integrate in [0, 1].\n           tmp = cis(x[1]*pi/2)*pi/2\n           # Assign to two components of \"f\" the real\n           # and imaginary part of the integrand.\n           f[1], f[2] = reim(tmp)\n       end\nintegrand (generic function with 1 method)\n\njulia> result = cuhre(integrand, 2, 2);\n\njulia> begin\n           println(\"Result of Cuba: \", complex(result[1]...))\n           println(\"Exact result:   \", complex(1.0, 1.0))\n       end\nResult of Cuba: 1.0 + 1.0im\nExact result:   1.0 + 1.0im"
},

{
    "location": "#Passing-data-to-the-integrand-function-1",
    "page": "Cuba",
    "title": "Passing data to the integrand function",
    "category": "section",
    "text": "Cuba Library allows program written in C and Fortran to pass extra data to the integrand function with userdata argument. This is useful, for example, when the integrand function depends on changing parameters. In Cuba.jl the userdata argument is not available, but you do not normally need it.For example, the cumulative distribution function F(xk) of chi-squared distribution is defined byF(x k) = int_0^x fract^k2 - 1exp(-t2)2^k2Gamma(k2)\nmathrmdtThe cumulative distribution function depends on parameter k, but the function passed as integrand to Cuba.jl integrator routines accepts as arguments only the input and output vectors. However you can easily define a function to calculate a numerical approximation of F(x k) based on the above integral expression because the integrand can access any variable visible in its scope. The following Julia script computes F(x = pi k) for different k and compares the result with more precise values, based on the analytic expression of the cumulative distribution function, provided by GSL.jl package.julia> using Cuba, GSL, Printf, SpecialFunctions\n\njulia> function chi2cdf(x::Real, k::Real)\n           k2 = k/2\n           # Chi-squared probability density function, without constant denominator.\n           # The result of integration will be divided by that factor.\n           function chi2pdf(t::Float64)\n               # \"k2\" is taken from the outside.\n               return t^(k2 - 1.0)*exp(-t/2)\n           end\n           # Neither \"x\" is passed directly to the integrand function,\n           # but is visible to it.  \"x\" is used to scale the function\n           # in order to actually integrate in [0, 1].\n           x*cuhre((t,f) -> f[1] = chi2pdf(t[1]*x))[1][1]/(2^k2*gamma(k2))\n       end\nchi2cdf (generic function with 1 method)\n\njulia> x = float(pi);\n\njulia> begin\n            @printf(\"Result of Cuba: %.6f %.6f %.6f %.6f %.6f\\n\",\n                    map((k) -> chi2cdf(x, k), collect(1:5))...)\n            @printf(\"Exact result:   %.6f %.6f %.6f %.6f %.6f\\n\",\n                    map((k) -> cdf_chisq_P(x, k), collect(1:5))...)\n        end\nResult of Cuba: 0.923681 0.792120 0.629694 0.465584 0.321833\nExact result:   0.923681 0.792120 0.629695 0.465584 0.321833"
},

{
    "location": "#Vectorized-Function-1",
    "page": "Cuba",
    "title": "Vectorized Function",
    "category": "section",
    "text": "Consider the integralintlimits_Omega prod_i=1^10 cos(x_i)\nmathrmdboldsymbolx = sin(1)^10 = 01779883dotswhere Omega = 0 1^10 and boldsymbolx = (x_1 dots x_10) is a 10-dimensional vector. A simple way to compute this integral is the following:julia> using Cuba, BenchmarkTools\n\njulia> cuhre((x, f) -> f[] = prod(cos.(x)), 10)\nComponent:\n 1: 0.1779870665870775 ± 1.0707995959536173e-6 (prob.: 0.2438374075714901)\nIntegrand evaluations: 7815\nFail:                  0\nNumber of subregions:  2\n\njulia> @benchmark cuhre((x, f) -> f[] = prod(cos.(x)), 10)\nBenchmarkTools.Trial:\n  memory estimate:  2.62 MiB\n  allocs estimate:  39082\n  --------------\n  minimum time:     1.633 ms (0.00% GC)\n  median time:      1.692 ms (0.00% GC)\n  mean time:        1.867 ms (8.62% GC)\n  maximum time:     3.660 ms (45.54% GC)\n  --------------\n  samples:          2674\n  evals/sample:     1We can use vectorization in order to speed up evaluation of the integrand function.julia> function fun_vec(x,f)\n           f[1,:] .= 1.0\n           for j in 1:size(x,2)\n               for i in 1:size(x, 1)\n                   f[1, j] *= cos(x[i, j])\n               end\n           end\n       end\nfun_vec (generic function with 1 method)\n\njulia> cuhre(fun_vec, 10, nvec = 1000)\nComponent:\n 1: 0.1779870665870775 ± 1.0707995959536173e-6 (prob.: 0.2438374075714901)\nIntegrand evaluations: 7815\nFail:                  0\nNumber of subregions:  2\n\njulia> @benchmark cuhre(fun_vec, 10, nvec = 1000)\nBenchmarkTools.Trial:\n  memory estimate:  2.88 KiB\n  allocs estimate:  54\n  --------------\n  minimum time:     949.976 μs (0.00% GC)\n  median time:      954.039 μs (0.00% GC)\n  mean time:        966.930 μs (0.00% GC)\n  maximum time:     1.204 ms (0.00% GC)\n  --------------\n  samples:          5160\n  evals/sample:     1A further speed up can be gained by running the for loop in parallel with Threads.@threads. For example, running Julia with 4 threads:julia> function fun_par(x,f)\n           f[1,:] .= 1.0\n           Threads.@threads for j in 1:size(x,2)\n               for i in 1:size(x, 1)\n                   f[1, j] *= cos(x[i, j])\n               end\n           end\n       end\nfun_par (generic function with 1 method)\n\njulia> cuhre(fun_par, 10, nvec = 1000)\nComponent:\n 1: 0.1779870665870775 ± 1.0707995959536173e-6 (prob.: 0.2438374075714901)\nIntegrand evaluations: 7815\nFail:                  0\nNumber of subregions:  2\n\njulia> @benchmark cuhre(fun_par, 10, nvec = 1000)\nBenchmarkTools.Trial:\n  memory estimate:  3.30 KiB\n  allocs estimate:  63\n  --------------\n  minimum time:     507.914 μs (0.00% GC)\n  median time:      515.182 μs (0.00% GC)\n  mean time:        520.667 μs (0.06% GC)\n  maximum time:     3.801 ms (85.06% GC)\n  --------------\n  samples:          9565\n  evals/sample:     1"
},

{
    "location": "#Performance-1",
    "page": "Cuba",
    "title": "Performance",
    "category": "section",
    "text": "Cuba.jl cannot (yet?) take advantage of parallelization capabilities of Cuba Library. Nonetheless, it has performances comparable with equivalent native C or Fortran codes based on Cuba library when CUBACORES environment variable is set to 0 (i.e., multithreading is disabled). The following is the result of running the benchmark present in test directory on a 64-bit GNU/Linux system running Julia 0.7.0-beta2.3 (commit 83ce9c7524) equipped with an Intel(R) Core(TM) i7-4700MQ CPU. The C and FORTRAN 77 benchmark codes have been compiled with GCC 7.3.0.$ CUBACORES=0 julia -e \'using Pkg; cd(Pkg.dir(\"Cuba\")); include(\"test/benchmark.jl\")\'\n[ Info: Performance of Cuba.jl:\n  0.257360 seconds (Vegas)\n  0.682703 seconds (Suave)\n  0.329552 seconds (Divonne)\n  0.233190 seconds (Cuhre)\n[ Info: Performance of Cuba Library in C:\n  0.268249 seconds (Vegas)\n  0.682682 seconds (Suave)\n  0.319553 seconds (Divonne)\n  0.234099 seconds (Cuhre)\n[ Info: Performance of Cuba Library in Fortran:\n  0.233532 seconds (Vegas)\n  0.669809 seconds (Suave)\n  0.284515 seconds (Divonne)\n  0.195740 seconds (Cuhre)Of course, native C and Fortran codes making use of Cuba Library outperform Cuba.jl when higher values of CUBACORES are used, for example:$ CUBACORES=1 julia -e \'using Pkg; cd(Pkg.dir(\"Cuba\")); include(\"test/benchmark.jl\")\'\n[ Info: Performance of Cuba.jl:\n  0.260080 seconds (Vegas)\n  0.677036 seconds (Suave)\n  0.342396 seconds (Divonne)\n  0.233280 seconds (Cuhre)\n[ Info: Performance of Cuba Library in C:\n  0.096388 seconds (Vegas)\n  0.574647 seconds (Suave)\n  0.150003 seconds (Divonne)\n  0.102817 seconds (Cuhre)\n[ Info: Performance of Cuba Library in Fortran:\n  0.094413 seconds (Vegas)\n  0.556084 seconds (Suave)\n  0.139606 seconds (Divonne)\n  0.107335 seconds (Cuhre)Cuba.jl internally fixes CUBACORES to 0 in order to prevent from forking julia processes that would only slow down calculations eating up the memory, without actually taking advantage of concurrency. Furthemore, without this measure, adding more Julia processes with addprocs() would only make the program segfault."
},

{
    "location": "#Related-projects-1",
    "page": "Cuba",
    "title": "Related projects",
    "category": "section",
    "text": "There are other Julia packages for multidimenensional numerical integration:Cubature.jl\nHCubature.jl\nNIntegration.jl"
},

{
    "location": "#Development-1",
    "page": "Cuba",
    "title": "Development",
    "category": "section",
    "text": "Cuba.jl is developed on GitHub: https://github.com/giordano/Cuba.jl. Feel free to report bugs and make suggestions at https://github.com/giordano/Cuba.jl/issues."
},

{
    "location": "#History-1",
    "page": "Cuba",
    "title": "History",
    "category": "section",
    "text": "The ChangeLog of the package is available in NEWS.md file in top directory. There have been some breaking changes from time to time, beware of them when upgrading the package."
},

{
    "location": "#License-1",
    "page": "Cuba",
    "title": "License",
    "category": "section",
    "text": "The Cuba.jl package is licensed under the GNU Lesser General Public License, the same as Cuba library. The original author is Mosè Giordano."
},

{
    "location": "#Credits-1",
    "page": "Cuba",
    "title": "Credits",
    "category": "section",
    "text": "If you use this library for your work, please credit Thomas Hahn. Citable papers about Cuba Library:Hahn, T. 2005, Computer Physics Communications, 168, 78. DOI:10.1016/j.cpc.2005.01.010. arXiv:hep-ph/0404043. Bibcode:2005CoPhC.168…78H.\nHahn, T. 2015, Journal of Physics Conference Series, 608, 012066. DOI:10.1088/1742-6596/608/1/012066. arXiv:1408.6373. Bibcode:2015JPhCS.608a2066H."
},

]}
