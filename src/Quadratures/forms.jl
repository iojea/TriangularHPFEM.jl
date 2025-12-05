"""
See `@form`

```
   struct Form
```

A struct for storing the list of integrands and measures that define a linear of bilinear form.

It is recommended to build `Form`s using the `@form` macro. Howvever, a `Form` can be created directly by passing a function to be integrated and a `Measure`.

```
    julia> Ω = circmesh(0.1);
    julia> dΩ = Measure(Ω);
    julia> Form((u,v)->∇(u)⋅∇(v),dΩ)
```

If the `Form` has two or more terms, tuples of functions to be integrated and its corresponding  `Measures` are expected.
```
   julia> term₁(u,v) = ∇(u)⋅∇(v);
   julia> term₂(u,v) = 2u*v;
   julia> Form((term₁,term₂),(dΩ,dΩ)) 
```
But when every term is integrated with respect to the same `Measure`, it can be passed once:
```
   julia> Form((term₁,term₂),dΩ) 
```
"""
struct Form
    nargs
    integrands
    measures
end

function Form(t::Tuple,m::Tuple)
    println(t)
    println(m)
    tnargs = Tuple(first(methods(tt)).nargs-1 for tt in t)
    println(tnargs)
    all(tnargs .== tnargs[1]) || throw(ArgumentError("Number of arguments vary from term to term."))
    all(typeof(mm)<:Measure for mm in m) || throw(ArgumentError("Measures are not measures."))
    println(typeof(t))
    println(typeof(m))
    Form(tnargs[1],t,m)
end
Form(f::Function,m::Measure) = Form((f,),(m,))
#Form(t::Tuple,m::Measure) = Form(t,Tuple(m for _ in 1:length(t)))

function _form(x...)
    n = length(x)
    n%2==0 || throw(ArgumentError("Malformed expression. Each term must be of the form `∫(fun)*dΩ` where `fun` is some function and `dΩ` is a `Measure`."))
    k = n÷2
    Form(tuple(x[1:k]...),tuple(x[k+1:end]...))
end
    
bad_integrand(::Any) = false
function bad_integrand(expr::Expr)
    expr.head == :call && expr.args[1] == :+ && return true
    expr.head == :call && expr.args[1] == :- && length(expr.args)>2 && return true
    expr.head == :call && return any(bad_integrand.(expr.args[2:end]))
    expr.head == :block && return any(bad_integrand.(expr.args))
    return false
end

function fracture(expr,arg,signed,others...)
    if expr.head == :call && expr.args[1] == :*
        measure = esc(expr.args[3])
        if expr.args[2].head == :call && expr.args[2].args[1] == :∫
            if signed 
                integrand = Expr(:call,:-,expr.args[2].args[2])
            else
                integrand = expr.args[2].args[2]
            end
            bad_integrand(integrand) && throw(error("Integrands containing sums are not supported. Please use one integral per term in your form. "))
            fun_body = Expr(:->,arg,integrand)
        elseif expr.args[2].head == :call && expr.args[2].args[1] in (:*,:⋅)
            throw(error("Multiplication outside of integrals are not supported yet. Place the multiplication inside of the integrand instead.")) 
        end
        return fun_body,measure,others...
    end
    if expr.head == :call && expr.args[1] == :+
        terms1 = fracture(expr.args[2],arg,false)
        terms2 = fracture(expr.args[3],arg,false)
        return terms1...,terms2...,others...
    elseif expr.head == :call && expr.args[1] == :-
        terms1 = fracture(expr.args[2],arg,false)
        terms2 = fracture(expr.args[3],arg,true)
        return terms1...,terms2...,others...
    end
end

"""
```
   @form form_definition
```
Creates a `Form` 
```
   julia> A = rand(2,2);
   julia> Ω = circmesh(0.1);
   julia> dΩ = Measure(Ω)
   julia> @form a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ
```
Creates a `Form` named `a`.

Linear `Forms` are built the same way
```
   julia> f(x) = x[1]^2*sin(x[2])
   julia> @form b(v) = ∫(v*f)*dΩ
```
It is also possible to create `Form`s with several terms.
```
    julia> @form a(u,v) = ∫(∇(u)*∇(v))*dΩ + ∫(3u*v)*dΩ
```
It is not allowed to have several terms under the same integral sign.
```
    julia> @form a(u,v)=∫(∇(u)⋅∇(v)+u*v)*dΩ
    ERROR: LoadError: Integrands containing sums are not supported. Please use one integral per term in your form.
```
This is mainly for performance reasons. The terms of first and second order are integrated differently, so they cannot be parts of the same function. 
""" 
macro form(expr)
       name = get_name(expr)
       par = get_parameters(expr)
       terms = get_terms(expr)
       tira = fracture(terms,par,false)
       Expr(:(=),name,Expr(:call,:_form,tira[1:2:end]...,tira[2:2:end]...))
end


get_name(expr) = esc(expr.args[1].args[1])
get_parameters(expr) = Expr(:tuple,expr.args[1].args[2:end]...)
get_terms(expr) = expr.args[2].args[2]


