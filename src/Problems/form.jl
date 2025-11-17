"""
```
   struct Form
```
A struct for storing the list of integrands and measures that define a linear of bilinear form.

It is recommended to build `Form`s using the `@form` macro:
```
   julia> A = rand(2,2);
   julia> Ω = circmesh(0.1);
   julia> dΩ = Measure(Ω)
   julia> @form a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ
```
Creates a `Form` named `a`

"""
struct Form
    nargs
    integrands
    measures
end
function Form(t::Tuple,m::Tuple) 
    tnargs = Tuple(methods(tt).ms[1].nargs-1 for tt in t)
    all(tnargs .== tnargs[1]) || throw(ArgumentError("Number of arguments vary from term to term."))
    all(typeof(mm)<:Measure for mm in m) || throw(ArgumentError("Measures are not measures."))
    Form(tnargs[1],t,m)
end
Form(f::Function,m::Measure) = Form((f,),(m,))


get_name(expr) = expr.args[1].args[1]
get_parameters(expr) = Expr(:tuple,expr.args[1].args[2:end]...)
get_integrand(expr) = expr.args[2].args[2]
get_terms(expr) = Meta.parse.(split(string(expr.args[2].args[2]),['+','-']))
get_measure(term) = term.args[3]
macro form(expr)
    par = get_parameters(expr)
    terms = get_terms(expr)
    codes = get_integrand.(terms)
    meas = get_measure.(terms)
    tira = Form(Tuple(eval(Expr(:->,par,code)) for code in codes),
                 Tuple(eval(m) for m in meas))
    name = get_name(expr)
    eval(Expr(:(=),name,:($tira)))
end



