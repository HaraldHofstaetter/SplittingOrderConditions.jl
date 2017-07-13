module SplittingOrderConditions

import Base: start, next, done

export Lyndon, MultiFor, multinomial_coeff, generate_equations, call_form
export generate_equations_ABC, get_ABC_coeffs
export compute_residua, compute_residua_ABC, get_lem, get_p_lem
   
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
                            last_b_is_zero::Bool=false, with_brackets::Bool=true)
    eqs = Dict{Array{Int,1}, AbstractString}()
    if symmetric && iseven(q)
        return eqs
    end
    if symmetric
        palindromic = false
        last_b_is_zero = true    
        s2 = ceil(Int,s/2)
        s3 = isodd(s)?s2-1:s2
        if with_brackets
            a = vcat([string("a[", j, "]") for j=1:s2], [string("a[", j, "]") for j=s3:-1:1])
        else    
            a = vcat([string("a", j) for j=1:s2], [string("a", j) for j=s3:-1:1])
        end
        s2 = ceil(Int,(s-1)/2)
        s3 = iseven(s)?s2-1:s2
        if with_brackets
            b = vcat([string("b[", j, "]") for j=1:s2], [string("b[", j, "]") for j=s3:-1:1],["0"])
        else
            b = vcat([string("b", j) for j=1:s2], [string("b", j) for j=s3:-1:1],["0"])
        end
    elseif palindromic
        last_b_is_zero = false 
        s2 = ceil(Int,s/2)
        s3 = isodd(s)?s2-1:s2
        if with_brackets
            a = vcat([string("a[", j, "]") for j=1:s2], [string("b[", j, "]") for j=s3:-1:1])
        else            
            a = vcat([string("a", j) for j=1:s2], [string("b", j) for j=s3:-1:1])
        end    
        b = reverse(a)
    else
        if with_brackets
            a = [string("a[", j, "]") for j=1:s]
            b = [string("b[", j, "]") for j=1:s]
        else    
            a = [string("a", j) for j=1:s]
            b = [string("b", j) for j=1:s]
        end
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


function generate_equations_ABC(q::Int, s::Int; symmetric::Bool=false, palindromic::Bool=false,
                            last_c_is_zero::Bool=false, last_b_is_zero::Bool=false, with_brackets::Bool=true)
    eqs = Dict{Array{Int,1}, AbstractString}()
    if symmetric && iseven(q)
        return eqs
    end
    if last_b_is_zero
       last_c_is_zero = true
    end   
    if symmetric
        error("symmertric not implemented.")
        #palindromic = false
        #last_c_is_zero = true    
        #last_b_is_zero = true    
        #s2 = ceil(Int,s/2)
        #s3 = isodd(s)?s2-1:s2
        #if with_brackets
        #    a = vcat([string("a[", j, "]") for j=1:s2], [string("a[", j, "]") for j=s3:-1:1])
        #else    
        #    a = vcat([string("a", j) for j=1:s2], [string("a", j) for j=s3:-1:1])
        #end
        #s2 = ceil(Int,(s-1)/2)
        #s3 = iseven(s)?s2-1:s2
        #if with_brackets
        #    b = vcat([string("b[", j, "]") for j=1:s2], [string("b[", j, "]") for j=s3:-1:1],["0"])
        #else
        #    b = vcat([string("b", j) for j=1:s2], [string("b", j) for j=s3:-1:1],["0"])
        #end
    elseif palindromic
        last_c_is_zero = false 
        last_b_is_zero = false 
        s2 = ceil(Int,s/2)
        s3 = isodd(s)?s2-1:s2
        if with_brackets
            a = vcat([string("a[", j, "]") for j=1:s2], [string("c[", j, "]") for j=s3:-1:1])
            b = vcat([string("b[", j, "]") for j=1:s2], [string("b[", j, "]") for j=s3:-1:1])
        else            
            a = vcat([string("a", j) for j=1:s2], [string("c", j) for j=s3:-1:1])
            b = vcat([string("b", j) for j=1:s2], [string("b", j) for j=s3:-1:1])
        end    
        c = reverse(a)
    else
        if with_brackets
            a = [string("a[", j, "]") for j=1:s]
            b = [string("b[", j, "]") for j=1:s]
            c = [string("c[", j, "]") for j=1:s]
        else    
            a = [string("a", j) for j=1:s]
            b = [string("b", j) for j=1:s]
            c = [string("c", j) for j=1:s]
        end
    end
    if palindromic
        rr = (2, 1, 0)
        for j in Lyndon(3, q)
            if length(j)==q && !haskey(eqs, [rr[x+1] for x in reverse(j)])
                eqs[j] = ""
            end
        end
    else
        for j in Lyndon(3, q)
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
        for lA = MultiFor(k)
        for lB = MultiFor(k-lA)
        for lC = MultiFor(k-lA-lB)
            W=Int[]
            for j=s:-1:1
                append!(W, 2*ones(Int, lC[j]))
                append!(W, ones(Int, lB[j]))
                append!(W, zeros(Int, lA[j]))
            end
            if W in keys(eqs) && !(last_b_is_zero && lB[s]>0)  && !(last_c_is_zero && lC[s]>0) 
                coeff =  multinomial_coeff(q, k)            
                for j=1:s      
                    coeff = coeff * multinomial_coeff(k[j], [lA[j], lB[j], lC[j]])
                end
                h = string(coeff)
                for j=1:s        
                    if lA[j]>0
                        h = string(h,"*",a[j])
                    end                    
                    if lA[j]>1
                        h = string(h,"^",lA[j])
                    end
                    if lB[j]>0
                        h = string(h,"*",b[j])
                    end    
                    if lB[j]>1
                        h = string(h,"^",lB[j])
                    end
                    if lC[j]>0
                        h = string(h,"*",c[j])
                    end    
                    if lC[j]>1
                        h = string(h,"^",lC[j])
                    end
                end
                eqs[W] = string(eqs[W],"+",h)
            end
        end
        end
        end
    end
    for (W, h) in eqs
        eqs[W] = string(h,"-1")
    end
    eqs
end


function compute_residua(q::Int, a::Vector, b::Vector)
    s = length(a)
    T = typeof(a[1])
    eqs = Dict{Array{Int,1}, T}()

    aa = zeros(T, s, q+1)
    bb = zeros(T, s, q+1)
    aa[:,1] = ones(T, s)
    bb[:,1] = ones(T, s)
    for k=1:q
        aa[:,k+1] = aa[:,k].*a
        bb[:,k+1] = bb[:,k].*b
    end

    for j in Lyndon(2, q)
        if length(j)==q
            eqs[j] = -one(T) 
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
            if W in keys(eqs)
                h =  multinomial_coeff(q, k)            
                for j=1:s      
                    h *= binomial(k[j], l[j])
                end
                for j=1:s        
                    h *= aa[j, k[j]-l[j]+1] * bb[j, l[j]+1] 
                end
                eqs[W] += h
            end
        end
    end
    eqs
end

get_lem(q, a, b) = norm(collect(values(compute_residua(q, a, b))))

function get_p_lem(a,b)
    q = 1    
    while true
        lem = get_lem(q,a,b) 
        if lem>1e-11
            return q-1, lem
        end
        q += 1
    end
end


function get_ABC_coeffs(a,b)
    s=length(a)
    T = typeof(a[1])
    S = 2*s-1
    A = zeros(T, S)
    B = zeros(T, S)
    C = zeros(T, S)
    h = a[1]
    for k=1:s
        K = 2*k-1
        A[K] = a[k]
        B[K] = h
        C[K] = b[k]
        h = b[k]-h
        if k<s
            A[K+1] = 0.0
            B[K+1] = h
            C[K+1] = 0.0
            h = a[k+1]-h
        end
    end
    A, B, C
end    



function compute_residua_ABC(q::Int, a::Vector, b::Vector, c::Vector)
    s = length(a)
    T = typeof(a[1])
    eqs = Dict{Array{Int,1}, T}()

    aa = zeros(T, s, q+1)
    bb = zeros(T, s, q+1)
    cc = zeros(T, s, q+1)
    aa[:,1] = ones(T, s)
    bb[:,1] = ones(T, s)
    cc[:,1] = ones(T, s)
    for k=1:q
        aa[:,k+1] = aa[:,k].*a
        bb[:,k+1] = bb[:,k].*b
        cc[:,k+1] = cc[:,k].*c
    end

    for j in Lyndon(3, q)
        if length(j)==q
            eqs[j] = -one(T) 
        end
    end

    for k0 in Combinatorics.WithReplacementCombinations(1:s,q)
        k = zeros(Int,s)
        for j in k0
            k[j]+=1
        end
        for lA = MultiFor(k)
        for lB = MultiFor(k-lA)
        for lC = MultiFor(k-lA-lB)
            W=Int[]
            for j=s:-1:1
                append!(W, 2*ones(Int, lC[j]))
                append!(W, ones(Int, lB[j]))
                append!(W, zeros(Int, lA[j]))
            end
            if W in keys(eqs)
                h =  multinomial_coeff(q, k)            
                for j=1:s      
                    h *= multinomial_coeff(k[j], [lA[j], lB[j], lC[j]])
                end
                for j=1:s        
                    h *= aa[j, lA[j]+1] * bb[j, lB[j]+1] * cc[j, lC[j]+1]
                end
                eqs[W] += h
            end
        end
        end
        end
    end
    eqs
end

get_lem(q, a, b, c) = norm(collect(values(compute_residua_ABC(q, a, b, c))))

function get_p_lem(a,b,c)
    q = 1    
    while true
        lem = get_lem(q,a,b,c) 
        if lem>1e-11
            return q-1, lem
        end
        q += 1
    end
end



function call_form(input::ASCIIString; threads=1)
    attempts = 0
    FORM_PID=0
    path = joinpath(dirname(@__FILE__), "../deps/bin")
    out = "Error initializing preset external channels\n"
    while out == "Error initializing preset external channels\n" 
        (so,si,pr) = readandwrite(`$path/tform -w$(threads) -t /tmp -M  -q -pipe 0,1 $path/pipe.frm`)
        FORM_PID = readuntil(so,'\n')[1:end-1];
        print(si, FORM_PID,',',getpid(),"\n\n", input,"\n\n")
        close(si)
        out = readall(so)
        close(si)
        close(so)
        close(pr)
        attempts +=1
    end
    println(STDERR, "# of attemps to call form = ",attempts)
    try 
        rm("/tmp/xform$FORM_PID.str") # delete generated temporary file
    end    
    out
end



function compile_fg(fun::AbstractString, n::Integer; threads=1)
    fun = replace(replace(fun,"[", "("), "]", ")")
    input = string("""
Off Statistics;
V x;
S u;
""","Local f0 = ",
    fun,
    ";\n",
"""
#write <> "begin"
#write <> "  if length(g)==0"
Print;
Format O3;
Format maple;
.sort;
#write <> "    return f0"
#write <> "  else"
I i;
S m, xx, u;
Hide f0;
#do i=1,$(n)
    Local f'i' = f0;
    id x('i') = xx;
    id xx^m? = m*xx^(m-1);
    id xx = x('i');
    .sort
    Hide f'i';
#enddo
.sort;
""","Local H = ",
    join(["u^$(j+1)*f$(j)" for j=0:n],:+),
    ";\n",
"""
B u;
.sort
#optimize H
B u;
.sort
""",
    join(["Local F$(j) = H[u^$(j+1)];\n" for j=0:n], ""),
"""
.sort
#write <> "%4O",
""", "#write <> \"\n   f=%e   ",
join(["g($(j))=%e" for j=1:n],"    "),
"\", ",
join(["F$(j)" for j=0:n],","),"\n",
"""
#write <> "    return f"
#write <> "  end"
#write <> "end"
.end
""")
    out = call_form(input, threads=threads)
    out = replace(replace(out,"(", "["), ")", "]")    
    out = replace(out,"length[g]", "length(g)")
    eval(parse(string("(x,g)->",out)))
end    



end
