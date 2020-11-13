"""
    module FMM3D

Wrappers for the Flatiron Institute's FMM3D
library. 

# wrappers

All N-body codes return output in an 
`FMMVals` structure. See documentation
of N-body codes for details.

N-body interactions with the Helmholtz kernel
- [`hfmm3d`](@ref): ``O(N)`` fast mutlipole code
- [`h3ddir`](@ref): ``O(N^2)`` direct code

N-body interactions with the Laplace kernel
- [`lfmm3d`](@ref): ``O(N)`` fast mutlipole code
- [`l3ddir`](@ref): ``O(N^2)`` direct code

# lower level routines

For a list of lower level routines see the 
documentation for [`lower_level_routs`](@ref)
"""
module FMM3D

using FMM3D_jll

export FMMVals, hfmm3d, lfmm3d, h3ddir, l3ddir, lowerlevel_routs
export besseljs3d

"""

# Available lower level routines

* [`besseljs3d`](@ref) spherical Bessel function eval
"""
function lower_level_routs() end

# fortran input/return types

Fd = Ref{Float64}
Fi = Ref{Int32}
Fc = Ref{ComplexF64}

# common input types

TFN = Union{Array{Float64},Nothing}
TCN = Union{Array{ComplexF64},Nothing}


"""
    FMMVals()

Return type for FMM and direct computation function
calls. Fields are `nothing` on return unless requested.
See individual FMM/direct computation function 
documentation for specifics.
"""
mutable struct FMMVals
    pot
    grad
    hess

    pottarg
    gradtarg
    hesstarg

    # stokes-specific return values
    
    pre
    pretarg

    # maxwell-specific return values
    e
    h
    grade
    gradh

    etarg
    htarg
    gradetarg
    gradhtarg
end

function FMMVals()
    return FMMVals(nothing,nothing,nothing,
                   nothing,nothing,nothing,
                   nothing,nothing,
                   nothing,nothing,nothing,nothing,
                   nothing,nothing,nothing,nothing)
end

"""
```julia
    vals = hfmm3d(eps,zk,sources;charges=nothing,dipvecs=nothing,
                  targets=nothing,pg=0,pgt=0,nd=1)
```

This function computes the N-body Helmholtz interactions
in three dimensions where the interaction kernel is given 
by ``e^{ikr}/r`` and its gradients. This is the 
``O(N)`` fast multipole code which computes the interactions
to the requested precision.

```math
 u(x) = \\sum_{j=1}^{N} c_{j} \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
- v_{j} \\cdot \\nabla \\left( \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
\\right)  \\, , 
```
 
where ``c_{j}`` are the charge densities,  
      ``v_{j}`` are the dipole orientation vectors, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum

# Input

* `eps::Float64` precision requested
* `zk::ComplexF64` Helmholtz parameter 
* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)
* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `pg::Integer` source eval flag. 
    + Potential (``u``) at sources evaluated if `pg == 1`. 
    + Potential and gradient (``\\nabla u``) at sources evaluated if `pg == 2`
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities

# Output
        
`vals::FMMVals` with the fields

* `vals.pot::Array{ComplexF64}` size (nd,n) or (n) potential at source locations if requested
* `vals.grad::Array{ComplexF64}` size (nd,3,n) or (3,n) gradient at source locations if requested
* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
"""
function hfmm3d(eps::Float64,zk::Union{Float64,ComplexF64},
                sources::Array{Float64};
                charges::TCN=nothing,dipvecs::TCN=nothing,
                targets::TFN=nothing,pg::Integer=0,pgt::Integer=0,
                nd::Integer=1)

    zk = complex(zk)
    
    # default values

    vals = FMMVals()
    
    ifcharge = 0
    ifdipole = 0

    zero = ComplexF64(0)
    
    pot = zero
    grad = zero
    hess = zero
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    
    @assert size(sources,1) == 3
    @assert nd >= 0

    if (pg > 2 || pg < 0); @warn "flag pg not in expected range" end
    if (pgt > 2 || pgt < 0); @warn "flag pgt not in expected range" end    

    n = div(length(sources),3)

    nt = 0

    if targets != nothing
        @assert size(targets,1) == 3
        nt = div(length(targets),3)
    else
        targets = 0.0
    end
    
    if charges != nothing
        @assert div(length(charges),nd) == n
        ifcharge = 1
    else
        charges = zero
        ifcharge = 0
    end

    if dipvecs != nothing
        @assert div(length(dipvecs),nd) == n*3
        ifdipole = 1
    else
        dipvecs = zero
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        return vals
    end

    if (pg != 1 && pg != 2) && (pgt != 1 && pgt != 2)
        @warn "no output requested, doing nothing"
        return vals
    end

    if (pgt == 1 || pgt == 2) && targets == nothing
        @warn "target values requested but no targets provided"
    end

    # allocate memory for return values
    
    if pg == 1 || pg == 2
        if nd > 1
            pot = zeros(ComplexF64,nd,n)
        else
            pot = zeros(ComplexF64,n)
        end
    end

    if pg == 2
        if nd > 1
            grad = zeros(ComplexF64,nd,3,n)
        else
            grad = zeros(ComplexF64,3,n)
        end
    end

    if pgt == 1 || pgt == 2
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,3,nt)
        else
            gradtarg = zeros(ComplexF64,3,nt)
        end
    end

    # actually call the function
    #
    # fortran calling sequence:
    # subroutine hfmm3d(nd,eps,zk,nsource,source,ifcharge,
    #     $    charge,ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,
    #     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg)
    #

    
    ccall((:hfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fc,Fi,Fd,Fi,Fc,
                                    Fi,Fc,Fi,Fc,Fc,Fc,Fi,
                                    Fd,Fi,Fc,Fc,Fc),
          nd,eps,zk,n,sources,ifcharge,charges,ifdipole,
          dipvecs,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg)

    # load requested values

    if pg == 1 || pg == 2; vals.pot = pot end
    if pg == 2; vals.grad = grad end
    if pgt == 1 || pgt == 2; vals.pottarg = pottarg end
    if pgt == 2; vals.gradtarg = gradtarg end
    
    return vals

end

"""
```julia
    vals = h3ddir(zk,sources,targets;charges=nothing,
                    dipvecs=nothing,pgt=0,nd=1,
                    thresh=1e-16)
```

This function computes the N-body Helmholtz interactions
in three dimensions where the interaction kernel is given 
by ``e^{ikr}/r`` and its gradients. This is the 
``O(N^2)`` direct evaluation code.

```math
 u(x) = \\sum_{j=1}^{N} c_{j} \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
- v_{j} \\cdot \\nabla \\left( \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
\\right)  \\, , 
```
 
where ``c_{j}`` are the charge densities,  
      ``v_{j}`` are the dipole orientation vectors, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum

# Input

* `zk::ComplexF64` Helmholtz parameter 
* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities
* `thresh::Float64` threshold for ignoring interactions when ``\\|x-x_{j}\\| \\leq thresh``

# Output
        
`vals::FMMVals` with the fields

* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
"""
function h3ddir(zk::Union{ComplexF64,Float64},sources::Array{Float64},
                targets::Array{Float64};
                charges::TCN=nothing,dipvecs::TCN=nothing,
                pgt::Integer=0,nd::Integer=1,
                thresh::Float64=1e-16)

    zk = complex(zk)
    
    # default values

    vals = FMMVals()
    
    ifcharge = 0
    ifdipole = 0

    zero = ComplexF64(0)
    
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    
    @assert size(sources,1) == 3
    @assert size(targets,1) == 3    
    @assert nd >= 0

    if (pgt > 2 || pgt < 0); @warn "flag pgt not in expected range" end    

    n = div(length(sources),3)
    nt = div(length(targets),3)
    
    if charges != nothing
        @assert div(length(charges),nd) == n
        ifcharge = 1
    else
        charges = zero
        ifcharge = 0
    end

    if dipvecs != nothing
        @assert div(length(dipvecs),nd) == n*3
        ifdipole = 1
    else
        dipvecs = zero
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        return vals
    end

    if (pgt != 1 && pgt != 2)
        @warn "no output requested, doing nothing"
        return vals
    end

    # allocate memory for return values
    
    if pgt == 1 || pgt == 2
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,3,nt)
        else
            gradtarg = zeros(ComplexF64,3,nt)
        end
    end

    # dispatch to appropriate wrapper
    
    if pgt == 1
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdp!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,thresh)
                vals.pottarg = pottarg
            else
                h3ddirectcp!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             thresh)
                vals.pottarg = pottarg
            end
        else
            if ifdipole == 1
                h3ddirectdp!(nd,zk,sources,
                             dipvecs,n,targets,nt,
                             pottarg,thresh)
                vals.pottarg = pottarg
            end
        end
    elseif pgt == 2
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdg!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                
            else
                h3ddirectcg!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                                
            end
        else
            if ifdipole == 1
                h3ddirectdg!(nd,zk,sources,
                              dipvecs,n,targets,nt,
                             pottarg,gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                                
            end
        end
    end                
            
    return vals

end

function h3ddirectcp!(nd,zk,sources,charges,n,targets,nt,
                      pot,thresh)
    ccall((:h3ddirectcp_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fd),
          nd,zk,sources,charges,n,targets,nt,pot,thresh)
    return
end

function h3ddirectdp!(nd,zk,sources,dipvecs,n,targets,nt,
                      pot,thresh)
    ccall((:h3ddirectdp_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fd),
          nd,zk,sources,dipvecs,n,targets,nt,pot,thresh)
    return
end

function h3ddirectcdp!(nd,zk,sources,charges,dipvecs,n,targets,
                       nt,pot,thresh)
    ccall((:h3ddirectcdp_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fc,Fi,
                                          Fd,Fi,Fc,Fd),
          nd,zk,sources,charges,dipvecs,n,targets,nt,pot,thresh)
    return
end

function h3ddirectcg!(nd,zk,sources,charges,n,targets,nt,
                      pot,grad,thresh)
    ccall((:h3ddirectcg_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fd),
          nd,zk,sources,charges,n,targets,nt,pot,grad,thresh)
    return
end

function h3ddirectdg!(nd,zk,sources,dipvecs,n,targets,nt,
                      pot,grad,thresh)
    ccall((:h3ddirectdg_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fd),
          nd,zk,sources,dipvecs,n,targets,nt,pot,grad,thresh)
    return
end

function h3ddirectcdg!(nd,zk,sources,charges,dipvecs,n,targets,
                       nt,pot,grad,thresh)
    ccall((:h3ddirectcdg_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fd),
          nd,zk,sources,charges,dipvecs,n,targets,nt,pot,grad,thresh)
    return
end

    
"""
```julia
    vals = lfmm3d(eps,sources;charges=nothing,dipvecs=nothing,
                  targets=nothing,pg=0,pgt=0,nd=1)
```

This function computes the N-body Laplace interactions
in three dimensions where the interaction kernel is given 
by ``\$1/r\$`` and its gradients.  This is the 
``O(N)`` fast multipole code which computes the interactions
to the requested precision.

```math
u(x) = \\sum_{j=1}^{N} c_{j} / \\|x-x_{j}\\| + v_{j} 
\\cdot \\nabla( 1/\\|x-x_{j}\\|)  \\, ,
u(x) = \\sum_{j=1}^{N} c_{j} \\frac{1}{\\|x-x_{j}\\|} 
- v_{j} \\cdot \\nabla \\left( \\frac{1}{\\|x-x_{j}\\|} 
\\right)  \\, , 
```
 
where ``c_{j}`` are the charge densities,  
      ``v_{j}`` are the dipole orientation vectors, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum

# Input

# Input

* `eps::Float64` precision requested
* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)
* `charges::Array{Float64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{Float64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `pg::Integer` source eval flag. 
    + Potential (``u``) at sources evaluated if `pg == 1`. 
    + Potential and gradient (``\\nabla u``) at sources evaluated if `pg == 2`
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities

# Output
        
`vals::FMMVals` with the fields

Note that the Hessian is returned as a length
6 vector at each point with the second derivatives
in the order: ``\\partial_{xx}``, ``\\partial_{yy}``,
``\\partial_{zz}``, ``\\partial_{xy}``, 
``\\partial_{xz}``, ``\\partial_{yz}``.

* `vals.pot::Array{Float64}` size (nd,n) potential at source locations if requested
* `vals.grad::Array{Float64}` size (nd,3,n) gradient at source locations if requested
* `vals.hess::Array{Float64}` size (nd,6,n) Hessian at source locations if requested
* `vals.pottarg::Array{Float64}` size (nd,nt) potential at target locations if requested
* `vals.gradtarg::Array{Float64}` size (nd,3,nt) gradient at target locations if requested
* `vals.hesstarg::Array{Float64}` size (nd,6,nt) Hessian at target locations if requested
"""
function lfmm3d(eps::Float64,sources::Array{Float64};charges::TFN=nothing,
                dipvecs::TFN=nothing,targets::TFN=nothing,
                pg::Integer=0,pgt::Integer=0,nd::Integer=1)

    # default values

    vals = FMMVals()
    
    ifcharge = 0
    ifdipole = 0

    zero = Float64(0)
    
    pot = zero
    grad = zero
    hess = zero
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    
    @assert size(sources,1) == 3
    @assert nd >= 0

    if (pg > 3 || pg < 0); @warn "flag pg not in expected range" end
    if (pgt > 3 || pgt < 0); @warn "flag pgt not in expected range" end    

    n = div(length(sources),3)

    nt = 0

    if targets != nothing
        @assert size(targets,1) == 3
        nt = div(length(targets),3)
    else
        targets = 0.0
    end
    
    if charges != nothing
        @assert div(length(charges),nd) == n
        ifcharge = 1
    else
        charges = zero
        ifcharge = 0
    end

    if dipvecs != nothing
        @assert div(length(dipvecs),nd) == n*3
        ifdipole = 1
    else
        dipvecs = zero
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        return vals
    end

    if ((pg != 1 && pg != 2 && pg !=3) &&
        (pgt != 1 && pgt != 2 && pgt != 3))
        @warn "no output requested, doing nothing"
        return vals
    end

    if (pgt == 1 || pgt == 2 || pgt == 3) && targets == nothing
        @warn "target values requested but no targets provided"
    end

    # allocate memory for return values
    
    if pg == 1 || pg == 2 || pg == 3
        if nd > 1
            pot = zeros(Float64,nd,n)
        else
            pot = zeros(Float64,n)
        end
    end

    if pg == 2 || pg == 3
        if nd > 1
            grad = zeros(Float64,nd,3,n)
        else
            grad = zeros(Float64,3,n)
        end
    end

    if pg == 3
        if nd > 1
            hess = zeros(Float64,nd,6,n)
        else
            hess = zeros(Float64,6,n)
        end
    end

    if pgt == 1 || pgt == 2 || pgt == 3
        if nd > 1
            pottarg = zeros(Float64,nd,nt)
        else
            pottarg = zeros(Float64,nt)
        end
    end

    if pgt == 2 || pgt == 3
        if nd > 1
            gradtarg = zeros(Float64,nd,3,nt)
        else
            gradtarg = zeros(Float64,3,nt)
        end
    end

    if pgt == 3
        if nd > 1
            hesstarg = zeros(Float64,nd,6,nt)
        else
            hesstarg = zeros(Float64,6,nt)
        end
    end

    # actually call the function
    #
    # fortran calling sequence:
    # subroutine lfmm3d(nd,eps,nsource,source,ifcharge,
    #     $    charge,ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,
    #     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg)
    #

    
    ccall((:lfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fi,Fd,Fi,Fd,
                                    Fi,Fd,Fi,Fd,Fd,Fd,Fi,
                                    Fd,Fi,Fd,Fd,Fd),
          nd,eps,n,sources,ifcharge,charges,ifdipole,
          dipvecs,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg)

    # load requested values

    if pg == 1 || pg == 2 || pg == 3; vals.pot = pot end
    if pg == 2 || pg == 3; vals.grad = grad end
    if pg == 3; vals.hess = hess end        
    if pgt == 1 || pgt == 2 || pgt == 3; vals.pottarg = pottarg end
    if pgt == 2 || pgt == 3; vals.gradtarg = gradtarg end
    if pgt == 3; vals.hesstarg = hesstarg end    
    
    return vals

end


"""
```julia
    vals = l3ddir(sources,targets;charges=nothing,
                dipvecs=nothing,pgt=0,nd=1,
                thresh=1e-16)
```

This function computes the N-body Laplace interactions
in three dimensions where the interaction kernel is given 
by ``1/r`` and its gradients.  This is the 
``O(N^2)`` direct evaluation code.

```math
 u(x) = \\sum_{j=1}^{N} c_{j} \\frac{1}{\\|x-x_{j}\\|} 
- v_{j} \\cdot \\nabla \\left( \\frac{1}{\\|x-x_{j}\\|} 
\\right)  \\, , 
```
 
where ``c_{j}`` are the charge densities,  
      ``v_{j}`` are the dipole orientation vectors, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum

# Input

* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities
* `thresh::Float64` threshold for ignoring interactions when ``\\|x-x_{j}\\| \\leq thresh``

# Output
        
`vals::FMMVals` with the fields

* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
"""
function l3ddir(sources::Array{Float64},targets::Array{Float64};
                charges::TFN=nothing,dipvecs::TFN=nothing,
                pgt::Integer=0,nd::Integer=1,
                thresh::Float64=1e-16)
    
    # default values

    vals = FMMVals()
    
    ifcharge = 0
    ifdipole = 0

    zero = Float64(0)
    
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    # check inputs
    
    @assert size(sources,1) == 3
    @assert size(targets,1) == 3    
    @assert nd >= 0

    if (pgt > 2 || pgt < 0); @warn "flag pgt not in expected range" end    

    n = div(length(sources),3)
    nt = div(length(targets),3)
    
    if charges != nothing
        @assert div(length(charges),nd) == n
        ifcharge = 1
    else
        charges = zero
        ifcharge = 0
    end

    if dipvecs != nothing
        @assert div(length(dipvecs),nd) == n*3
        ifdipole = 1
    else
        dipvecs = zero
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        return vals
    end

    if (pgt != 1 && pgt != 2)
        @warn "no output requested, doing nothing"
        return vals
    end

    # allocate memory for return values
    
    if pgt == 1 || pgt == 2
        if nd > 1
            pottarg = zeros(Float64,nd,nt)
        else
            pottarg = zeros(Float64,nt)
        end
    end

    if pgt == 2
        if nd > 1
            gradtarg = zeros(Float64,nd,3,nt)
        else
            gradtarg = zeros(Float64,3,nt)
        end
    end

    # dispatch to appropriate wrapper
    
    if pgt == 1
        if ifcharge == 1
            if ifdipole == 1
                l3ddirectcdp!(nd,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,thresh)
                vals.pottarg = pottarg
            else
                l3ddirectcp!(nd,sources,charges,
                             n,targets,nt,pottarg,
                             thresh)
                vals.pottarg = pottarg
            end
        else
            if ifdipole == 1
                l3ddirectdp!(nd,sources,
                             dipvecs,n,targets,nt,
                             pottarg,thresh)
                vals.pottarg = pottarg
            end
        end
    elseif pgt == 2
        if ifcharge == 1
            if ifdipole == 1
                l3ddirectcdg!(nd,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                
            else
                l3ddirectcg!(nd,sources,charges,
                             n,targets,nt,pottarg,
                             gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                                
            end
        else
            if ifdipole == 1
                l3ddirectdg!(nd,sources,
                              dipvecs,n,targets,nt,
                             pottarg,gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                                
            end
        end
    end                
            
    return vals

end

function l3ddirectcp!(nd,sources,charges,n,targets,nt,
                      pot,thresh)
    ccall((:l3ddirectcp_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fi,
                                          Fd,Fi,Fd,Fd),
          nd,sources,charges,n,targets,nt,pot,thresh)
    return
end

function l3ddirectdp!(nd,sources,dipvecs,n,targets,nt,
                      pot,thresh)
    ccall((:l3ddirectdp_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fi,
                                          Fd,Fi,Fd,Fd),
          nd,sources,dipvecs,n,targets,nt,pot,thresh)
    return
end

function l3ddirectcdp!(nd,sources,charges,dipvecs,n,targets,
                       nt,pot,thresh)
    ccall((:l3ddirectcdp_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fd,Fi,
                                          Fd,Fi,Fd,Fd),
          nd,sources,charges,dipvecs,n,targets,nt,pot,thresh)
    return
end

function l3ddirectcg!(nd,sources,charges,n,targets,nt,
                      pot,grad,thresh)
    ccall((:l3ddirectcg_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fi,
                                          Fd,Fi,Fd,Fd,Fd),
          nd,sources,charges,n,targets,nt,pot,grad,thresh)
    return
end

function l3ddirectdg!(nd,sources,dipvecs,n,targets,nt,
                      pot,grad,thresh)
    ccall((:l3ddirectdg_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fi,
                                          Fd,Fi,Fd,Fd,Fd),
          nd,sources,dipvecs,n,targets,nt,pot,grad,thresh)
    return
end

function l3ddirectcdg!(nd,sources,charges,dipvecs,n,targets,
                       nt,pot,grad,thresh)
    ccall((:l3ddirectcdg_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fd,Fi,
                                          Fd,Fi,Fd,Fd,Fd),
          nd,sources,charges,dipvecs,n,targets,nt,pot,grad,thresh)
    return
end

"""
    function fj, fjder = besseljs3d(nterms,z;scale=1.0,ifder=0)

This subroutine evaluates the first `nterms` spherical Bessel 
functions, and if requested, their derivatives.
It incorporates a scaling parameter `scale` so that
      
      	fjs_n(z)=j_n(z)/SCALE^n
      	fjder_n(z)=\\frac{\\partial fjs_n(z)}{\\partial z}

# Input

* `nterms::Integer` order of expansion of output array `fjs` 
* `z::ComplexF64` argument of the spherical Bessel functions
* `scale::Float64` scaling factor
* `ifder::Integer1` flag indicating whether to calculate `fjder`
      	          0	NO
      	          1	YES
OUTPUT:

* `fjs::Array{ComplexF64}` array of length `nterms+1` of scaled Bessel functions.
* `fjder::Array{ComplexF64}` array of derivatives of scaled Bessel functions, if requested.
"""
function besseljs3d(nterms::Integer,z::ComplexF64;scale::Float64=1.0,
                    ifder::Integer=0)

    @assert (nterms >= 0)
    if (ifder !=0 && ifder != 1)
        @warn "unexpected value in ifder, no ders computed"
    end

    fjs = Array{ComplexF64}(undef,nterms+1)
    fjder = ComplexF64(0)

    
    
    if ifder == 1
        fjder = Array{ComplexF64}(undef,nterms+1)
    end
    
    # fortran interface
    # subroutine besseljs3d(nterms,z,scale,fjs,ifder,fjder)
    
    ccall((:besseljs3d_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,Fc),
          nterms,z,scale,fjs,ifder,fjder)

    if ifder != 1
        fjder = nothing
    end

    return fjs, fjder
    
end

end # module
