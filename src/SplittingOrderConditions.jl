module SplittingOrderConditions

import Base: start, next, done

export Lyndon, MultiFor, multinomial_coeff, generate_equations
   
using Combinatorics


immutable Lyndon
    s::Int
    n::Int
end

Base.start(::Lyndon) = Int[-1]

function Base.next(L::Lyndon, w::Vector{Int})
    if w == [-1]
        w[end] += 1
        return (copy(w), w)
    end
    m = length(w)
    while length(w) < L.n               
        push!(w, w[end-m+1])
    end
    while length(w) > 0 && w[end] == L.s - 1    
        pop!(w)
    end
    w[end] += 1
    (copy(w), w)
end
    
Base.done(L::Lyndon, w::Vector{Int}) = w == [L.s-1]

immutable MultiFor
    k::Array{Int,1}
end

Base.start(MF::MultiFor) = Int[]

Base.done(MF::MultiFor, k::Array{Int,1}) = MF.k==k

function Base.next(MF::MultiFor, k::Array{Int,1}) 
    if k==Int[]
        k = zeros(Int, length(MF.k))
        return(copy(k), k)
    end
    for i=1:length(k)
        if k[i]<MF.k[i]
            k[i] += 1
            for j = 1:i-1                 
                k[j] = 0       
            end
            return (copy(k), k)
        end
    end            
end


multinomial_coeff(q::Int, k::Array{Int,1}) = div(factorial(q), prod([factorial(i) for i in k]))



function generate_equations(q::Int, s::Int; symmetric::Bool=false, palindromic::Bool=false,
                            last_b_is_zero::Bool=false)
    eqs = Dict{Array{Int,1}, AbstractString}()
    if symmetric && iseven(q)
        return eqs
    end
    if symmetric
        palindromic = false
        last_b_is_zero = true    
        s2 = ceil(Int,s/2)
        s3 = isodd(s)?s2-1:s2
        a = vcat([string("a[", j, "]") for j=1:s2], [string("a[", j, "]") for j=s3:-1:1])
        s2 = ceil(Int,(s-1)/2)
        s3 = iseven(s)?s2-1:s2
        b = vcat([string("b[", j, "]") for j=1:s2], [string("b[", j, "]") for j=s3:-1:1],["0"])
    elseif palindromic
        last_b_is_zero = false 
        s2 = ceil(Int,s/2)
        s3 = isodd(s)?s2-1:s2
        a = vcat([string("a[", j, "]") for j=1:s2], [string("b[", j, "]") for j=s3:-1:1])
        b = reverse(a)
    else
        a = [string("a[", j, "]") for j=1:s]
        b = [string("b[", j, "]") for j=1:s]
    end
    if palindromic
        for j in Lyndon(2, q)
            if length(j)==q && !haskey(eqs, [x==0?1:0 for x in reverse(j)])
                eqs[j] = ""
            end
        end
    else
        for j in Lyndon(2, q)
            if length(j)==q
                eqs[j] = ""
            end
        end
    end    

    for k0 in Combinatorics.WithReplacementCombinations(1:s,q)
        k = zeros(Int,s)
        for j in k0
            k[j]+=1
        end
        for l = MultiFor(k)
            W=Int[]
            for j=s:-1:1
                append!(W, ones(Int, l[j]))
                append!(W, zeros(Int, k[j]-l[j]))
            end
            if W in keys(eqs) && !(last_b_is_zero && l[s]>0)
                coeff =  multinomial_coeff(q, k)            
                for j=1:s      
                    coeff = coeff * binomial(k[j], l[j])
                end
                h = string(coeff)
                for j=1:s        
                    if k[j]-l[j]>0
                        h = string(h,"*",a[j])
                    end                    
                    if k[j]-l[j]>1
                        h = string(h,"^",k[j]-l[j])
                    end
                    if l[j]>0
                        h = string(h,"*",b[j])
                    end    
                    if l[j]>1
                        h = string(h,"^",l[j])
                    end
                end
                eqs[W] = string(eqs[W],"+",h)
            end
        end
    end
    for (W, h) in eqs
        eqs[W] = string(h,"-1")
    end
    eqs
end

end
