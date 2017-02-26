module SplittingOrderConditions

import Base: start, next, done

export Lyndon, MultiFor, multinomial_coeff, generate_equations, call_form
   
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
join(["F$(j)" for j=0:n],","),
"""
#write <> "    return f"
#write <> "  end"
#write <> "end"
.end
""")
    out = call_form(input, threads=threads)
    out = replace(replace(out,"(", "["), ")", "]")    
    eval(parse(string("(x,g)->",out)))
end    



end
