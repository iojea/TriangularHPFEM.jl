function Base.show(io::IO,mime::MIME"text/plain",m::Measure) 
    println("Measure")
    header = Markdown.parse("""
        + Mesh with $(length(m.mesh.points)) nodes.
        + Scheme of degree $(m.sch.degree).
        """)
    show(io,mime,header)
end

function Base.show(io::IO,mime::MIME"text/plain",f::Form)
    if f.nargs == 2
        text = Markdown.parse("Bivariate form")
        show(io,mime,text)
    elseif f.nargs == 1
        text = Markdown.parse("Univariate form")
        show(io,mime,text)
    end
end