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
            return (copy(k), k)
        end
    end            
end


multinomial_coeff(q::Int, k::Array{Int,1}) = div(factorial(q), prod([factorial(i) for i in k]))



function generate_equations(q::Int, s::Int)
    L = Dict{Array{Int,1}, AbstractString}()
    for i in Lyndon(2, q)
        if length(i)==q
            L[i] = ""
        end
    end

    for k0 in Combinatorics.WithReplacementCombinations(1:s,q)
        k = zeros(Int,s)
        for j in k0
            k[j]+=1
        end
        for l = MultiFor(k)
            K=Int[]
            for j=s:-1:1
                append!(K, ones(Int, l[j]))
                append!(K, zeros(Int, k[j]-l[j]))
            end
            if K in keys(L)
                coeff =  multinomial_coeff(q, k)            
                for j=1:s      
                    coeff = coeff * binomial(k[j], l[j])
                end
                h = string(coeff)
                for j=1:s            
                    if k[j]-l[j]>0
                        h = string(h,"*a[",j,"]")
                    end                    
                    if k[j]-l[j]>1
                        h = string(h,"^",k[j]-l[j])
                    end
                    if l[j]>0
                        h = string(h,"*b[",j,"]")
                    end    
                    if l[j]>1
                        h = string(h,"^",l[j])
                    end
                end
                L[K] = string(L[K],"+",h)
            end
        end
    end
    for (K, h) in L
        L[K] = string(h,"-1=0")
    end
    L
end

end