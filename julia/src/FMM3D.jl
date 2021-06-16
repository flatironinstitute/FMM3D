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

N-body interactions with Stokes kernels
- [`stfmm3d`](@ref): ``O(N)`` fast mutlipole code
- [`st3ddir`](@ref): ``O(N^2)`` direct code

N-body interactions with an electromagnetic kernel
- [`emfmm3d`](@ref): ``O(N)`` fast mutlipole code
- [`em3ddir`](@ref): ``O(N^2)`` direct code

# lower level routines

For a list of lower level routines see the 
documentation for [`lower_level_routs`](@ref)
"""
module FMM3D

using FMM3D_jll

export FMMVals, hfmm3d, lfmm3d, h3ddir, l3ddir, lowerlevel_routs
export stfmm3d, st3ddir, emfmm3d, em3ddir
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

    ier

    # stokes-specific return values
    
    pre
    pretarg

    # maxwell-specific return values
    E
    divE
    curlE

    Etarg
    divEtarg
    curlEtarg
end

function FMMVals()
    return FMMVals(nothing,nothing,nothing,
                   nothing,nothing,nothing,
                   nothing,
                   nothing,nothing,
                   nothing,nothing,nothing,
                   nothing,nothing,nothing)
end


"""
```julia
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm3dinputcheck(sources,charges,dipvecs,
                                targets,pg,pgt,nd) )
```

Input checking for Helmholtz and Laplace routines. 
Checks the sizes of these arrays and the compatibility
of the various flags, nd, and provided arrays. 
If something is off, a warning is issued and 
`anyfail` is set to true. Non-fatal mistakes result in
a warning but `anyfail` remains false. 

Output:
* `anyfail` - boolean is true if any checks fail, false otherwise
* n - integer, number of sources
* nt - integer, number of targets
* ifcharge - integer, is 1 if there are charges, 0 otherwise
* ifdipole - integer, is 1 if there are dipoles, 0 otherwise
"""
function scalarfmm3dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd)
    anyfail = false
    
    if (size(sources,1) != 3)
        @warn "sources array is wrong shape, computing nothing"
        anyfail = true
    end
    if (nd < 0)
        @warn "nd has invalid value, computing nothing"
        anyfail = true
    end
    if (pg > 3 || pg < 0)
        @warn "flag pg not in expected range, computing nothing"
        anyfail = true
    end
    if (pgt > 3 || pgt < 0)
        @warn "flag pgt not in expected range, computing nothing"
        anyfail = true
    end    

    n = div(length(sources),3)
    nt = 0

    if targets != nothing
        if (size(targets,1) != 3)
            @warn "targets array is wrong shape, computing nothing"
            anyfail = true
        end
        nt = div(length(targets),3)
    end
    
    if charges != nothing
        if (div(length(charges),nd) != n)
            @warn "size of charges array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifcharge = 1
    else
        ifcharge = 0
    end

    if dipvecs != nothing
        if (div(length(dipvecs),nd) != n*3)
            @warn "size of dipvecs array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifdipole = 1
    else
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        anyfail = true
    end

    if (pg != 1 && pg != 2 && pg != 3) && (pgt != 1 && pgt != 2 && pgt != 3)
        @warn "no output requested, doing nothing"
        anyfail = true
    end

    if (pgt == 1 || pgt == 2 || pgt == 3) && targets == nothing
        @warn "target values requested but no targets provided, proceeding anyway"
    end

    return anyfail, n, nt, ifcharge, ifdipole
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

# Optional Keyword Input

* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `pg::Integer` source eval flag. 
    + Do not compute any quantities at sources if `pg == 0`
    + Potential (``u``) at sources evaluated if `pg == 1`. 
    + Potential and gradient (``\\nabla u``) at sources evaluated if `pg == 2`
    + Potential, gradient, and Hessian (``\\nabla \\nabla u``) at sources evaluated if `pg == 3`
* `pgt::Integer` target eval flag. 
    + Do not compute any quantities at targets if `pgt == 0`
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
    + Potential, gradient, and Hessian at targets evaluated if `pgt == 3`
* `nd::Integer` number of densities

Note: if all default values are used for optional input, nothing
is computed.

# Output
        
`vals::FMMVals` with the fields

Note that the Hessian is returned as a length
6 vector at each point with the second derivatives
in the order: ``\\partial_{xx}``, ``\\partial_{yy}``,
``\\partial_{zz}``, ``\\partial_{xy}``, 
``\\partial_{xz}``, ``\\partial_{yz}``.

* `vals.pot::Array{ComplexF64}` size (nd,n) or (n) potential at source locations if requested
* `vals.grad::Array{ComplexF64}` size (nd,3,n) or (3,n) gradient at source locations if requested
* `vals.hess::Array{ComplexF64}` size (nd,6,n) or (6,n) Hessian at source locations if requested
* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
* `vals.hesstarg::Array{ComplexF64}` size (nd,6,nt) or (6,nt) Hessian at target locations if requested
* `vals.ier` error flag as returned by FMM3D library. A value of 0 indicates a successful call. 
Non-zero values may indicate insufficient memory available. See the documentation for the FMM3D library. 
If not set (`nothing`), then FMM3D library was never called.
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

    if (pg == 3)
        @warn "Hessian not implemented for Helmholtz fmm, only computing potential and gradients at sources"
        pg = 2
    end
    if (pgt == 3)
        @warn "Hessian not implemented for Helmholtz fmm, only computing potential and gradients at targets"
        pgt = 2
    end
    
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm3dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if (ifcharge == 0); charges = zero end
    if (ifdipole == 0); dipvecs = zero end
    if (nt == 0); targets = 0.0 end
    
    # allocate memory for return values
    
    if pg == 1 || pg == 2 || pg == 3
        if nd > 1
            pot = zeros(ComplexF64,nd,n)
        else
            pot = zeros(ComplexF64,n)
        end
    end

    if pg == 2 || pg == 3
        if nd > 1
            grad = zeros(ComplexF64,nd,3,n)
        else
            grad = zeros(ComplexF64,3,n)
        end
    end

    if pg == 3
        if nd > 1
            hess = zeros(ComplexF64,nd,6,n)
        else
            hess = zeros(ComplexF64,6,n)
        end
    end

    if pgt == 1 || pgt == 2 || pgt == 3
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2 || pgt == 3
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,3,nt)
        else
            gradtarg = zeros(ComplexF64,3,nt)
        end
    end

    if pgt == 3
        if nd > 1
            hesstarg = zeros(ComplexF64,nd,6,nt)
        else
            hesstarg = zeros(ComplexF64,6,nt)
        end
    end

    # actually call the function
    #
    # fortran calling sequence:
    # subroutine hfmm3d(nd,eps,zk,nsource,source,ifcharge,
    #     $    charge,ifdipole,dipvec,iper,ifpgh,pot,grad,hess,ntarg,
    #     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
    #

    ier = Integer(0)
    iper = Integer(0)
    
    ccall((:hfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fc,Fi,Fd,Fi,Fc,
                                    Fi,Fc,Fi,Fi,Fc,Fc,Fc,Fi,
                                    Fd,Fi,Fc,Fc,Fc,Fi),
          nd,eps,zk,n,sources,ifcharge,charges,ifdipole,
          dipvecs,iper,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg,ier)

    vals.ier = ier
    if (ier != 0); @warn "libfmm3d had an error, see vals.ier"; end    
    
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
    vals = h3ddir(zk,sources,targets;charges=nothing,
                    dipvecs=nothing,pgt=0,nd=1,
                    thresh=1e-16)
```

This function computes the N-body Helmholtz interactions
in three dimensions where the interaction kernel is given 
by ``e^{ikr}/r`` and its gradients. This is the 
``O(N^2)`` direct evaluation code. By convention this code only computes 
the effect of sources on targets. If the value at sources is 
also needed, the routine can be called again with targets equal
to the source locations.

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

# Optional Keyword Input

* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
    + Potential, gradient, and Hessian at targets evaluated if `pgt == 3`
* `nd::Integer` number of densities
* `thresh::Float64` threshold for ignoring interactions when ``\\|x-x_{j}\\| \\leq thresh``

Note: if all default values are used for optional input, nothing
is computed.

# Output
        
`vals::FMMVals` with the fields

Note that the Hessian is returned as a length
6 vector at each point with the second derivatives
in the order: ``\\partial_{xx}``, ``\\partial_{yy}``,
``\\partial_{zz}``, ``\\partial_{xy}``, 
``\\partial_{xz}``, ``\\partial_{yz}``.

* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
* `vals.hesstarg::Array{ComplexF64}` size (nd,6,nt) or (6,nt) gradient at target locations if requested
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

    if (pgt == 3)
        @warn "Hessian not implemented for Helmholtz fmm, only computing potential and gradients at targets"
        pgt = 2
    end

    
    pg = 0
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm3dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if (ifcharge == 0); charges = zero end
    if (ifdipole == 0); dipvecs = zero end
    if (nt == 0); targets = 0.0 end

    # allocate memory for return values
    
    if pgt == 1 || pgt == 2 || pgt == 3
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2 || pgt == 3
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,3,nt)
        else
            gradtarg = zeros(ComplexF64,3,nt)
        end
    end

    if pgt == 3
        if nd > 1
            hesstarg = zeros(ComplexF64,nd,6,nt)
        else
            hesstarg = zeros(ComplexF64,6,nt)
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
    elseif pgt == 3
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdh!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,gradtarg,hesstarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg
                vals.hesstarg = hesstarg                                
            else
                h3ddirectch!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             gradtarg,hesstarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg
                vals.hesstarg = hesstarg                                                
            end
        else
            if ifdipole == 1
                h3ddirectdh!(nd,zk,sources,
                              dipvecs,n,targets,nt,
                             pottarg,gradtarg,hesstarg,
                             thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg
                vals.hesstarg = hesstarg                                                
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

function h3ddirectch!(nd,zk,sources,charges,n,targets,nt,
                      pot,grad,hess,thresh)
    ccall((:h3ddirectch_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fc,Fd),
          nd,zk,sources,charges,n,targets,nt,pot,grad,hess,thresh)
    return
end

function h3ddirectdh!(nd,zk,sources,dipvecs,n,targets,nt,
                      pot,grad,hess,thresh)
    ccall((:h3ddirectdh_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fc,Fd),
          nd,zk,sources,dipvecs,n,targets,nt,pot,grad,hess,thresh)
    return
end

function h3ddirectcdh!(nd,zk,sources,charges,dipvecs,n,targets,
                       nt,pot,grad,hess,thresh)
    ccall((:h3ddirectcdh_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fc,Fd),
          nd,zk,sources,charges,dipvecs,n,targets,nt,pot,grad,hess,thresh)
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

* `eps::Float64` precision requested
* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)

# Optional Keyword Input

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

Note: if all default values are used for optional input, nothing
is computed.

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
* `vals.ier` error flag as returned by FMM3D library. A value of 0 indicates a successful call. 
Non-zero values may indicate insufficient memory available. See the documentation for the FMM3D library. 
If not set (`nothing`), then FMM3D library was never called.
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

    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm3dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if (ifcharge == 0); charges = zero end
    if (ifdipole == 0); dipvecs = zero end
    if (nt == 0); targets = 0.0 end

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
    #     $    charge,ifdipole,dipvec,iper,ifpgh,pot,grad,
    #     $    hess,ntarg,
    #     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg,ier)
    #


    ier = Integer(0)
    iper = Integer(0)
    
    ccall((:lfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fi,Fd,Fi,Fd,
                                    Fi,Fd,Fi,Fi,Fd,Fd,Fd,Fi,
                                    Fd,Fi,Fd,Fd,Fd,Fi),
          nd,eps,n,sources,ifcharge,charges,ifdipole,
          dipvecs,iper,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg,ier)

    if (ier != 0); @warn "libfmm3d had an error, see vals.ier"; end
    vals.ier = ier

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
``O(N^2)`` direct evaluation code. By convention this code only computes 
the effect of sources on targets. If the value at sources is 
also needed, the routine can be called again with targets equal
to the source locations.

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

# Optional Keyword Input

* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities
* `thresh::Float64` threshold for ignoring interactions when ``\\|x-x_{j}\\| \\leq thresh``

Note: if all default values are used for optional input, nothing
is computed.

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


    pg = 0
    anyfail, n, nt, ifcharge, ifdipole = (
        scalarfmm3dinputcheck(sources,charges,dipvecs,targets,pg,pgt,nd))

    if anyfail
        return vals
    end

    if (ifcharge == 0); charges = zero end
    if (ifdipole == 0); dipvecs = zero end
    if (nt == 0); targets = 0.0 end

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
            else
                l3ddirectcp!(nd,sources,charges,
                             n,targets,nt,pottarg,
                             thresh)
            end
        else
            if ifdipole == 1
                l3ddirectdp!(nd,sources,
                             dipvecs,n,targets,nt,
                             pottarg,thresh)
            end
        end
    elseif pgt == 2
        if ifcharge == 1
            if ifdipole == 1
                l3ddirectcdg!(nd,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,gradtarg,thresh)
            else
                l3ddirectcg!(nd,sources,charges,
                             n,targets,nt,pottarg,
                             gradtarg,thresh)
            end
        else
            if ifdipole == 1
                l3ddirectdg!(nd,sources,
                              dipvecs,n,targets,nt,
                             pottarg,gradtarg,thresh)
            end
        end
    end                

    if pgt == 1 || pgt == 2 ; vals.pottarg = pottarg end
    if pgt == 2; vals.gradtarg = gradtarg end
    
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

    if (nterms < 0)
        @warn "nterms has invalid value, computing nothing"
        fjs = nothing; fjder = nothing;
        return fjs, fjder
    end
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


"""
```julia
    anyfail, n, nt, ifstoklet, ifstrslet = (
        stfmm3dinputcheck(sources,stoklet,strslet,strsvec,
                                targets,ppreg,ppregt,nd) )
```

Input checking for Stokes routines. 
Checks the sizes of these arrays and the compatibility
of the various flags, nd, and provided arrays. 
If something is off, a warning is issued and 
`anyfail` is set to true. Non-fatal mistakes result in
a warning but `anyfail` remains false. 

Output:
* `anyfail` - boolean is true if any checks fail, false otherwise
* n - integer, number of sources
* nt - integer, number of targets
* ifstoklet - integer, is 1 if there are Stokeslets, 0 otherwise
* ifstrslet - integer, is 1 if there are (type I) stresslets, 0 otherwise
"""
function stfmm3dinputcheck(sources,stoklet,strslet,strsvec,targets,ppreg,ppregt,nd)
    anyfail = false
    
    if (size(sources,1) != 3)
        @warn "sources array is wrong shape, computing nothing"
        anyfail = true
    end
    if (nd < 0)
        @warn "nd has invalid value, computing nothing"
        anyfail = true
    end
    if (ppreg > 3 || ppreg < 0)
        @warn "flag ppreg not in expected range, computing nothing"
        anyfail = true
    end
    if (ppregt > 3 || ppregt < 0)
        @warn "flag ppregt not in expected range, computing nothing"
        anyfail = true
    end    

    n = div(length(sources),3)
    nt = 0

    if targets != nothing
        if (size(targets,1) != 3)
            @warn "targets array is wrong shape, computing nothing"
            anyfail = true
        end
        nt = div(length(targets),3)
    end
    
    if stoklet != nothing
        if (div(length(stoklet),nd) != 3*n)
            @warn "size of stoklet array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifstoklet = 1
    else
        ifstoklet = 0
    end

    if (strslet != nothing || strsvec != nothing)
        if (strsvec == nothing || strslet == nothing)
            @warn "if stresslets in calculation, need both strslet and strsvec arrays, computing nothing"
            anyfail = true
        end
        if (div(length(strslet),nd) != n*3)
            @warn "size of strslet array is incompatible with sources array and nd paramter, computing nothing"
            anyfail = true
        end
        if (div(length(strsvec),nd) != n*3)
            @warn "size of strsvec array is incompatible with sources array and nd paramter, computing nothing"
            anyfail = true
        end
        ifstrslet = 1
    else
        ifstrslet = 0
    end

    if ifstoklet == 0 && ifstrslet == 0
        @warn "no Stokeslets or stresslets provided, doing nothing"
        anyfail = true
    end

    if (ppreg != 1 && ppreg != 2 && ppreg != 3) && (ppregt != 1 && ppregt != 2 && ppregt != 3)
        @warn "no output requested, doing nothing"
        anyfail = true
    end

    if (ppregt == 1 || ppregt == 2 || ppregt == 3) && targets == nothing
        if !anyfail
            @warn "target values requested but no targets provided, proceeding anyway"
        else
            @warn "target values requested but no targets provided. There are other points of failure. Computing nothing."
        end
    end

    return anyfail, n, nt, ifstoklet, ifstrslet
end


"""
```julia
    vals = stfmm3d(eps,sources;stoklet=nothing,strslet=nothing,
                   strsvec=nothing,targets=nothing,ppreg=0,
                   ppregt=0,nd=1)
```

This function computes the N-body Stokes interactions
in three dimensions where the interaction kernels are
the Stokeslet and stresslet (see below). This is the 
``O(N)`` fast multipole code which computes the interactions
to the requested precision.

We take the following conventions for the Stokes kernels:

For a source ``y`` and target ``x``, let ``r_i = x_i-y_i``
 and let ``r = \\sqrt{r_1^2 + r_2^2 + r_3^2}``

The Stokeslet, ``G_{ij}``, and its associated pressure tensor, 
``P_j``, (without the ``1/4\\pi`` scaling) are

```math
G_{ij}(x,y) = (r_i r_j)/(2r^3) + \\delta_{ij}/(2r) \\; ,
\\quad
P_j(x,y) = r_j/r^3
```
The (Type I) stresslet, ``T_{ijk}``, and its associated 
pressure tensor, ``\\Pi_{jk}``, (without the ``1/4\\pi``
 scaling) are
         
```math
 T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5 \\; , \\quad
 PI_{jk} = -2 \\delta_{jk} + 6 r_j r_k/r^5   
```

The output of this routine gives the velocity
```math
 u_i(x) = \\sum_{m=1}^{N} G_{ij}(x,x_m) \\sigma^{(m)}_j 
    + T_{ijk}(x,x_m) \\mu^{(m)}_j \\nu^{(m)}_k \\, , 
```
and, optionally, the pressure

```math
 p(x) = \\sum_{m=1}^{N} P_{j}(x,x_m) \\sigma^{(m)}_j 
    + \\Pi{jk}(x,x_m) \\mu^{(m)}_j \\nu^{(m)}_k \\, , 
```
where ``\\sigma^{(m)}`` is a Stokeslet density 3-vector,  
``\\mu^{(m)}`` and ``\\nu^{(m)}`` are the stresslet 
density and orientation 3-vectors, and
``x_{m}`` are the source locations.
When ``x=x_{m}``, the term corresponding to 
``x_{m}`` is dropped from the sum.
The gradient of the velocity can also be computed by request.

# Input

* `eps::Float64` precision requested
* `sources::Array{Float64}` size (3,n) source locations (``x_{m}``)

# Optional Keyword Input

* `stoklet::Array{Float64}` size (nd,3,n) or (3,n) Stokeslet densities (``\\sigma``)
* `strslet::Array{Float64}` size (nd,3,n) or (3,n) Stresslet densities (``\\mu'')  if either strsvec or strslet is provided, the other must be as well
* `strsvec::Array{Float64}` size (nd,3,n) or (3,n) Stresslet orientations (``\\nu``) if either strsvec or strslet is provided, the other must be as well
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `ppreg::Integer` source eval flag. 
    + Velocity (``u``) at sources evaluated if `ppreg == 1`. 
    + Velocity and pressure (``p``) at sources evaluated if `ppreg == 2`
    + Velocity, pressure, and velocity gradient (``\\nabla u``) 
    at sources evaluated if `ppreg == 3`
* `ppregt::Integer` source eval flag. 
    + Velocity (``u``) at targets evaluated if `ppregt == 1`. 
    + Velocity and pressure (``p``) at targets evaluated if `ppregt == 2`
    + Velocity, pressure, and velocity gradient (``\\nabla u``) 
    at targets evaluated if `ppregt == 3`
* `nd::Integer` number of densities of each type

Note: if all default values are used for optional input, nothing
is computed.

# Output
        
`vals::FMMVals` with the fields

* `vals.pot::Array{Float64}` size (nd,3,n) or (3,n) velocity at source locations if requested
* `vals.grad::Array{Float64}` size (nd,3,3,n) or (3,n) gradient of velocity at source locations if requested.
`grad[l,i,j,k]` is the i-th component of the gradient of the j-th component of the velocity
for the l-th density at the k-th source location
* `vals.pre::Array{Float64}` size (nd,n) or (n) pressure at source locations if requested
* `vals.pottarg::Array{Float64}` size (nd,3,nt) or (3,nt) velocity at target locations if requested
* `vals.gradtarg::Array{Float64}` size (nd,3,3,nt) or (3,nt) gradient of velocity at target locations if requested.
`gradtarg[l,i,j,k]` is the i-th component of the gradient of the j-th component of the velocity
for the l-th density at the k-th target location
* `vals.pretarg::Array{Float64}` size (nd,nt) or (nt) pressure at target locations if requested
* `vals.ier` error flag as returned by FMM3D library. A value of 0 indicates a successful call. 
Non-zero values may indicate insufficient memory available. See the documentation for the FMM3D library. 
If not set (`nothing`), then FMM3D library was never called.
"""
function stfmm3d(eps::Float64,
                sources::Array{Float64};
                 stoklet::TFN=nothing,strslet::TFN=nothing,
                 strsvec::TFN=nothing,
                 targets::TFN=nothing,ppreg::Integer=0,
                 ppregt::Integer=0,nd::Integer=1)

    # default values

    vals = FMMVals()
    
    ifstoklet = 0
    ifstrslet = 0

    zero = Float64(0)
    zero3 = zeros(Float64,3)
    
    pot = zero
    grad = zero
    pre = zero
    pottarg = zero
    gradtarg = zero
    pretarg = zero

    ier = Integer(0)

    # check inputs

    anyfail, n, nt, ifstoklet, ifstrslet = (
        stfmm3dinputcheck(sources,stoklet,strslet,strsvec,targets,ppreg,ppregt,nd))

    if anyfail
        return vals
    end

    if (ifstoklet == 0); stoklet = zero3 end
    if (ifstrslet == 0); strslet = zero3; strsvec=zero3 end
    if (nt == 0); targets = zero3 end

    # allocate memory for return values
    
    if ppreg == 1 || ppreg == 2 || ppreg == 3
        if nd > 1
            pot = zeros(Float64,nd,3,n)
        else
            pot = zeros(Float64,3,n)
        end
    end

    if ppreg == 2 || ppreg == 3
        if nd > 1
            pre = zeros(Float64,nd,n)
        else
            pre = zeros(Float64,n)
        end
    end

    if ppreg == 3
        if nd > 1
            grad = zeros(Float64,nd,3,3,n)
        else
            grad = zeros(Float64,3,3,n)
        end
    end

    if ppregt == 1 || ppregt == 2 || ppregt == 3
        if nd > 1
            pottarg = zeros(Float64,nd,3,nt)
        else
            pottarg = zeros(Float64,3,nt)
        end
    end

    if ppregt == 2 || ppregt == 3
        if nd > 1
            pretarg = zeros(Float64,nd,nt)
        else
            pretarg = zeros(Float64,nt)
        end
    end

    if ppregt == 3
        if nd > 1
            gradtarg = zeros(Float64,nd,3,3,nt)
        else
            gradtarg = zeros(Float64,3,3,nt)
        end
    end

    # actually call the function
    #
    # fortran calling sequence
    #  subroutine stfmm3d(nd, eps, 
    # $                 nsource, source,
    # $                 ifstoklet, stoklet, ifstrslet, strslet, strsvec,
    # $                 ifppreg, pot, pre, grad, ntarg, targ, 
    # $                 ifppregtarg, pottarg, pretarg, gradtarg,ier)    

    
    ccall((:stfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fi,Fd,Fi,Fd,
                                    Fi,Fd,Fd,Fi,Fd,Fd,Fd,Fi,
                                    Fd,Fi,Fd,Fd,Fd,Fi),
          nd,eps,n,sources,ifstoklet,stoklet,ifstrslet,
          strslet,strsvec,ppreg,pot,pre,grad,nt,targets,ppregt,
          pottarg,pretarg,gradtarg,ier)

    # copy over error message

    if (ier != 0); @warn "libfmm3d had an error, see vals.ier" end
    vals.ier = ier
    
    # load requested values

    if ppreg == 1 || ppreg == 2 || ppreg == 3; vals.pot = pot end
    if ppreg == 2 || ppreg == 3; vals.pre = pre end
    if ppreg == 3; vals.grad = grad end    
    if ppregt == 1 || ppregt == 2 || ppregt == 3; vals.pottarg = pottarg end
    if ppregt == 2 || ppregt == 3; vals.pretarg = pretarg end
    if ppregt == 3; vals.gradtarg = gradtarg end    

    
    return vals

end


"""
```julia
    vals = st3ddir(sources,targets;stoklet=nothing,strslet=nothing,
                   strsvec=nothing,ppregt=0,nd=1,thresh=1e-16)
```

This function computes the N-body Stokes interactions
in three dimensions where the interaction kernels are
the Stokeslet and stresslet (see below). This is the 
``O(N^2)`` direct evaluation code. By convention this code only computes 
the effect of sources on targets. If the value at sources is 
also needed, the routine can be called again with targets equal
to the source locations.

We take the following conventions for the Stokes kernels:

For a source ``y`` and target ``x``, let ``r_i = x_i-y_i``
 and let ``r = \\sqrt{r_1^2 + r_2^2 + r_3^2}``

The Stokeslet, ``G_{ij}``, and its associated pressure tensor, 
``P_j``, (without the ``1/4\\pi`` scaling) are

```math
G_{ij}(x,y) = (r_i r_j)/(2r^3) + \\delta_{ij}/(2r) \\; ,
\\quad
P_j(x,y) = r_j/r^3
```
The (Type I) stresslet, ``T_{ijk}``, and its associated 
pressure tensor, ``\\Pi_{jk}``, (without the ``1/4\\pi``
 scaling) are
         
```math
 T_{ijk}(x,y) = -3 r_i r_j r_k/ r^5 \\; , \\quad
 PI_{jk} = -2 \\delta_{jk} + 6 r_j r_k/r^5   
```

The output of this routine gives the velocity
```math
 u_i(x) = \\sum_{m=1}^{N} G_{ij}(x,x_m) \\sigma^{(m)}_j 
    + T_{ijk}(x,x_m) \\mu^{(m)}_j \\nu^{(m)}_k \\, , 
```
and, optionally, the pressure

```math
 p(x) = \\sum_{m=1}^{N} P_{j}(x,x_m) \\sigma^{(m)}_j 
    + \\Pi{jk}(x,x_m) \\mu^{(m)}_j \\nu^{(m)}_k \\, , 
```
where ``\\sigma^{(m)}`` is a Stokeslet density 3-vector,  
``\\mu^{(m)}`` and ``\\nu^{(m)}`` are the stresslet 
density and orientation 3-vectors, and
``x_{m}`` are the source locations.
When ``x=x_{m}``, the term corresponding to 
``x_{m}`` is dropped from the sum.
The gradient of the velocity can also be computed by request.

# Input

* `sources::Array{Float64}` size (3,n) source locations (``x_{m}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)

# Optional Keyword Input

* `stoklet::Array{Float64}` size (nd,3,n) or (3,n) Stokeslet densities (``\\sigma``)
* `strslet::Array{Float64}` size (nd,3,n) or (3,n) Stresslet densities (``\\mu'')  if either strsvec or strslet is provided, the other must be as well
* `strsvec::Array{Float64}` size (nd,3,n) or (3,n) Stresslet orientations (``\\nu``) if either strsvec or strslet is provided, the other must be as well
* `ppregt::Integer` source eval flag. 
    + Velocity (``u``) at targets evaluated if `ppregt == 1`. 
    + Velocity and pressure (``p``) at targets evaluated if `ppregt == 2`
    + Velocity, pressure, and velocity gradient (``\\nabla u``) 
    at targets evaluated if `ppregt == 3`
* `nd::Integer` number of densities of each type

Note: if all default values are used for optional input, nothing
is computed.

# Output
        
`vals::FMMVals` with the fields

* `vals.pottarg::Array{Float64}` size (nd,3,nt) or (3,nt) velocity at target locations if requested
* `vals.gradtarg::Array{Float64}` size (nd,3,3,nt) or (3,nt) gradient of velocity at target locations if requested.
`gradtarg[l,i,j,k]` is the i-th component of the gradient of the j-th component of the velocity
for the l-th density at the k-th target location
* `vals.pretarg::Array{Float64}` size (nd,nt) or (nt) pressure at target locations if requested
"""
function st3ddir(sources::Array{Float64},targets::Array{Float64};
                 stoklet::TFN=nothing,strslet::TFN=nothing,
                 strsvec::TFN=nothing,
                 ppregt::Integer=0,nd::Integer=1,thresh=1e-16)

    # default values

    vals = FMMVals()
    
    ifstoklet = 0
    ifstrslet = 0

    zero = Float64(0)
    zero3 = zeros(Float64,3)
    
    pottarg = zero
    gradtarg = zero
    pretarg = zero

    ier = Integer(0)

    # check inputs

    ppreg=0
    anyfail, n, nt, ifstoklet, ifstrslet = (
        stfmm3dinputcheck(sources,stoklet,strslet,strsvec,targets,ppreg,ppregt,nd))

    if anyfail
        return vals
    end

    if (ifstoklet == 0); stoklet = zero3 end
    if (ifstrslet == 0); strslet = zero3 end
    if (nt == 0); targets = zero3 end

    # allocate memory for return values

    # for now we always allocate all quantities
    # and return requested
    
    if nd > 1
        pottarg = zeros(Float64,nd,3,nt)
    else
        pottarg = zeros(Float64,3,nt)
    end

    if nd > 1
        pretarg = zeros(Float64,nd,nt)
    else
        pretarg = zeros(Float64,nt)
    end

    if nd > 1
        gradtarg = zeros(Float64,nd,3,3,nt)
    else
        gradtarg = zeros(Float64,3,3,nt)
    end

    # dispatch to appropriate wrapper
    
    if ifstoklet == 1
        if ifstrslet == 1
            istress = 1
            st3ddirectstokstrsg!(nd,sources,stoklet,
                                 istress,strslet,strsvec,
                                 n,targets,nt,pottarg,
                                 pretarg,gradtarg,thresh)
        else
            st3ddirectstokg!(nd,sources,stoklet,
                             n,targets,nt,pottarg,
                             pretarg,gradtarg,thresh)
        end
    else
        if ifstrslet == 1
            istress = 1
            stoklet = zeros(Float64,nd,3,n)
            st3ddirectstokstrsg!(nd,sources,stoklet,
                             istress,strslet,strsvec,
                             n,targets,nt,pottarg,
                             pretarg,gradtarg,thresh)
        end
    end
    
    if ppregt == 1 || ppregt == 2 || ppregt == 3; vals.pottarg = pottarg end
    if ppregt == 2 || ppregt == 3; vals.pretarg = pretarg end
    if ppregt == 3; vals.gradtarg = gradtarg end    

    
    return vals

end


function st3ddirectstokstrsg!(nd,sources,stoklet,istress,
                              strslet,strsvec,n,targets,nt,pottarg,
                              pretarg,gradtarg,thresh)
    ccall((:st3ddirectstokstrsg_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fi,Fd,Fd,Fi,
                                                  Fd,Fi,Fd,Fd,Fd,Fd),
          nd,sources,stoklet,istress,
                              strslet,strsvec,n,targets,nt,pottarg,
          pretarg,gradtarg,thresh)
    return
end

function st3ddirectstokg!(nd,sources,stoklet,n,targets,nt,pottarg,
                              pretarg,gradtarg,thresh)
    ccall((:st3ddirectstokg_,libfmm3d),Cvoid,(Fi,Fd,Fd,Fi,Fd,Fi,Fd,
                                              Fd,Fd,Fd),
          nd,sources,stoklet,n,targets,nt,pottarg,
          pretarg,gradtarg,thresh)
    return
end

"""
```julia
    anyfail, n, nt, ifA, ifB, iflam = (
        emfmm3dinputcheck(sources,A,B,lambda,
                        targets,ifE,ifdivE,ifcurlE,
                        ifEtarg,ifdivEtarg,
                        ifcurlEtarg,nd) )
```

Input checking for electromagnetics routines. 
Checks the sizes of these arrays and the compatibility
of the various flags, nd, and provided arrays. 
If something is off, a warning is issued and 
`anyfail` is set to true. Non-fatal mistakes result in
a warning but `anyfail` remains false. 

Output:
* `anyfail` - boolean is true if any checks fail, false otherwise
* `n` - integer, number of sources
* `nt` - integer, number of targets
* `ifA` - integer, is 1 if there is an A vector, 0 otherwise
* `ifB` - integer, is 1 if there is a B vector, 0 otherwise
* `iflam` - integer, is 1 if there is a lambda density, 0 otherwise
"""
function emfmm3dinputcheck(sources,A,B,lambda,targets,
                           ifE,ifdivE,ifcurlE,ifEtarg,ifdivEtarg,
                           ifcurlEtarg,nd)
    anyfail = false
    
    if (size(sources,1) != 3)
        @warn "sources array is wrong shape, computing nothing"
        anyfail = true
    end
    if (nd < 0)
        @warn "nd has invalid value, computing nothing"
        anyfail = true
    end

    n = div(length(sources),3)
    nt = 0

    if targets != nothing
        if (size(targets,1) != 3)
            @warn "targets array is wrong shape, computing nothing"
            anyfail = true
        end
        nt = div(length(targets),3)
    end
    
    if A != nothing
        if (div(length(A),nd) != 3*n)
            @warn "size of A vector array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifA = 1
    else
        ifA = 0
    end

    if B != nothing
        if (div(length(B),nd) != 3*n)
            @warn "size of B vector array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        ifB = 1
    else
        ifB = 0
    end

    if lambda != nothing
        if (div(length(lambda),nd) != n)
            @warn "size of lambda array incompatible with sources array and nd parameter, computing nothing"
            anyfail = true
        end
        iflam = 1
    else
        iflam = 0
    end

    if ifA == 0 && ifB == 0 && iflam == 0
        @warn "none of A vectors, B vectors, or lambdas provided at sources, doing nothing"
        anyfail = true
    end

    if !ifEtarg && !ifdivEtarg && !ifcurlEtarg
        @warn "no output requested, doing nothing"
        anyfail = true
    end

    if (ifEtarg || ifdivEtarg || ifcurlEtarg) && targets == nothing
        @warn "target values requested but no targets provided, doing nothing for these"
    end

    return anyfail, n, nt, ifA, ifB, iflam
end


"""
```julia
    vals = emfmm3d(eps,zk,sources;A=nothing,B=nothing,lambda=nothing,
                ifE=false,ifdivE=false,ifcurlE=false,
                ifEtarg=false,ifdivEtarg=false,ifcurlEtarg=false,
                nd=1,targets=nothing)
```

This function computes N-body electromagnetic interactions
in three dimensions where the interaction kernels are
the Helmholtz kernel and its curl and gradient (see below). This is the 
``O(N)`` fast multipole code which computes the interactions
to the requested precision.

The Helmholtz Green's function (without the ``1/4 \\pi`` scaling)
is
```math
 G_k(x,y) = \\frac{e^{ik \\|x-y\\|}}{\\|x-y\\|} \\, , 
```

This routine computes the sum 

```math
E(x) = \\sum_{j=1}^m \\nabla \\times (G_k(x,x_{j}) A_j) + 
G_k(x,x_{j}) B_j + \\nabla G_k(x,x_{j}) \\lambda_j
```
 
where ``A_{j}`` and ``B_{j}`` are vector densities
      ``\\lambda_{j}`` is a scalar density, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum


# Input

* `eps::Float64` precision requested
* `zk::ComplexF64` Helmholtz parameter
* `sources::Array{Float64}` size (3,n) source locations (``x_{m}``)

# Optional Keyword Input

* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `A::Array{ComplexF64}` size (nd,3,n) or (3,n) vector densities
* `B::Array{ComplexF64}` size (nd,3,n) or (3,n) vector densities
* `lambda::Array{ComplexF64}` size (nd,n) or (n) scalar densities
* `ifE::Bool` E field evaluation flag, if true evaluate E at source points. 
* `ifdivE::Bool` divergence of E field evaluation flag, if true evaluate divergence of E at source points. 
* `ifcurlE::Bool` curl of E field evaluation flag, if true evaluate curl of E at source points.
* `ifEtarg::Bool` E field evaluation flag, if true evaluate E at target points. 
* `ifdivEtarg::Bool` divergence of E field evaluation flag, if true evaluate divergence of E at target points. 
* `ifcurlEtarg::Bool` curl of E field evaluation flag, if true evaluate curl of E at target points. 
* `nd::Integer` number of densities of each type

Note: if all default values are used for optional input, nothing
is computed.

# Output
        
`vals::FMMVals` with the fields

* `vals.E::Array{ComplexF64}` size (nd,3,n) or (3,n) E field at source locations if requested
* `vals.divE::Array{ComplexF64}` size (nd,n) or (n) divergence of E field at source locations if requested
* `vals.curlE::Array{ComplexF64}` size (nd,3,n) or (3,n) curl of E field at source locations if requested
* `vals.Etarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) E field at target locations if requested
* `vals.divEtarg::Array{ComplexF64}` size (nd,nt) or (nt) divergence of E field at target locations if requested
* `vals.curlEtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) curl of E field at target locations if requested
* `vals.ier` error flag as returned by FMM3D library. A value of 0 indicates a successful call. 
Non-zero values may indicate insufficient memory available. See the documentation for the FMM3D library. 
If not set (`nothing`), then FMM3D library was never called.
"""
function emfmm3d(eps::Float64,zk::Union{Float64,ComplexF64},sources::Array{Float64};
                 A::TCN=nothing,B::TCN=nothing,lambda::TCN=nothing,
                 ifE::Bool=false,ifdivE::Bool=false,
                 ifcurlE::Bool=false,ifEtarg::Bool=false,
                 ifdivEtarg::Bool=false,ifcurlEtarg::Bool=false,
                 nd::Integer=1,targets::TFN=nothing)

    zk = complex(zk)
    
    # default values

    vals = FMMVals()
    
    zero = ComplexF64(0)
    zero3 = zeros(ComplexF64,3)
    
    E1 = zero3
    divE1 = zero
    curlE1 = zero3

    ier = Integer(0)

    # check inputs

    anyfail, n, nt, ifA, ifB, iflam = (
        emfmm3dinputcheck(sources,A,B,lambda,
                          targets,ifE,ifdivE,ifcurlE,
                          ifEtarg,ifdivEtarg,
                          ifcurlEtarg,nd) )

    if anyfail
        return vals
    end

    if (ifA == 0); A = zero3 end
    if (ifB == 0); B = zero3 end
    if (iflam == 0); lambda = zero end    
    if (nt == 0); targets = zero3 end

    # allocate memory for return values

    anysrc = (ifE || ifdivE || ifcurlE)
    anytarg = (ifEtarg || ifdivEtarg || ifcurlEtarg)

    ntot = 0
    if anysrc; ntot = ntot+n; end
    if anytarg; ntot = ntot+nt; end

    targ1 = zeros(Float64,3,ntot)
    itstart = 1
    itend = nt
    if anysrc
        targ1[:,1:n] = sources
        itstart = itstart+n
        itend = itend+n
    end
    if anytarg
        targ1[:,itstart:itend] = targets
    end

    anyE = ifE || ifEtarg
    anydivE = ifdivE || ifdivEtarg
    anycurlE = ifcurlE || ifcurlEtarg

    ifE1 = 0; ifdivE1 = 0; ifcurlE1 = 0
    
    if anyE
        if nd > 1
            E1 = zeros(ComplexF64,nd,3,ntot)
        else
            E1 = zeros(ComplexF64,3,ntot)
        end
        ifE1 = 1
    end

    if anydivE
        if nd > 1
            divE1 = zeros(ComplexF64,nd,ntot)
        else
            divE1 = zeros(ComplexF64,ntot)
        end
        ifdivE1 = 1
    end

    if anyE
        if nd > 1
            curlE1 = zeros(ComplexF64,nd,3,ntot)
        else
            curlE1 = zeros(ComplexF64,3,ntot)
        end
        ifcurlE1 = 1
    end


# actually call the function
#
# fortran calling sequence
#  subroutine emfmm3d(nd,eps,zk,ns,source,ifa_vect,a_vect,&
#   ifb_vect,b_vect,iflambda,lambda,nt,targets,ifE,E,ifcurlE,curlE,&
#   ifdivE,divE,ier)

    
    ccall((:emfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fc,Fi,Fd,Fi,Fc,
                                      Fi,Fc,Fi,Fc,Fi,Fd,Fi,Fc,
                                      Fi,Fc,Fi,Fc,Fi),
          nd,eps,zk,n,sources,ifA,A,ifB,B,iflam,lambda,
          ntot,targ1,ifE1,E1,ifcurlE1,curlE1,ifdivE1,divE1,ier)

    # copy over error message

    if (ier != 0); @warn "libfmm3d had an error, see vals.ier" end
    vals.ier = ier
    
    # load requested values

    if (nd == 1)
        if ifE; vals.E = E1[:,1:n] end
        if ifdivE; vals.divE = divE1[1:n] end
        if ifcurlE; vals.curlE = curlE1[:,1:n] end
        
        if ifEtarg; vals.Etarg = E1[:,itstart:itend] end
        if ifdivEtarg; vals.divEtarg = divE1[itstart:itend] end
        if ifcurlEtarg; vals.curlEtarg = curlE1[:,itstart:itend] end
    else
        if ifE; vals.E = E1[:,:,1:n] end
        if ifdivE; vals.divE = divE1[:,1:n] end
        if ifcurlE; vals.curlE = curlE1[:,:,1:n] end
        
        if ifEtarg; vals.Etarg = E1[:,:,itstart:itend] end
        if ifdivEtarg; vals.divEtarg = divE1[:,itstart:itend] end
        if ifcurlEtarg; vals.curlEtarg = curlE1[:,:,itstart:itend] end
    end    
    return vals

end


"""
```julia
    vals = em3ddir(zk,sources,targets;A=nothing,B=nothing,lambda=nothing,
                ifEtarg=false,ifdivEtarg=false,ifcurlEtarg=false,
                nd=1,thresh=1e-16)
```

This function computes N-body electromagnetic interactions
in three dimensions where the interaction kernels are
the Helmholtz kernel and its curl and gradient (see below). This is the 
``O(N^2)`` direct code. By convention this code only computes 
the effect of sources on targets. If the value at sources is 
also needed, the routine can be called again with targets equal
to the source locations.

The Helmholtz Green's function (without the ``1/4 \\pi`` scaling)
is
```math
 G_k(x,y) = \\frac{e^{ik \\|x-y\\|}}{\\|x-y\\|} \\, , 
```

This routine computes the sum 

```math
E(x) = \\sum_{j=1}^m \\nabla \\times (G_k(x,x_{j}) A_j) + 
G_k(x,x_{j}) B_j + \\nabla G_k(x,x_{j}) \\lambda_j
```
 
where ``A_{j}`` and ``B_{j}`` are vector densities
      ``\\lambda_{j}`` is a scalar density, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum


# Input

* `zk::ComplexF64` Helmholtz parameter
* `sources::Array{Float64}` size (3,n) source locations (``x_{m}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)

# Optional Keyword Input

* `A::Array{ComplexF64}` size (nd,3,n) or (3,n) vector densities
* `B::Array{ComplexF64}` size (nd,3,n) or (3,n) vector densities
* `lambda::Array{ComplexF64}` size (nd,n) or (n) scalar densities
* `ifEtarg::Bool` E field evaluation flag, if true evaluate E at target points. 
* `ifdivEtarg::Bool` divergence of E field evaluation flag, if true evaluate divergence of E at target points. 
* `ifcurlEtarg::Bool` curl of E field evaluation flag, if true evaluate curl of E at target points. 
* `nd::Integer` number of densities of each type

Note: if all default values are used for optional input, nothing
is computed.

# Output
        
`vals::FMMVals` with the fields

* `vals.Etarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) E field at target locations if requested
* `vals.divEtarg::Array{ComplexF64}` size (nd,nt) or (nt) divergence of E field at target locations if requested
* `vals.curlEtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) curl of E field at target locations if requested
"""
function em3ddir(zk::Union{Float64,ComplexF64},sources::Array{Float64},
                 targets::Array{Float64};
                 A::TCN=nothing,B::TCN=nothing,lambda::TCN=nothing,
                 ifEtarg::Bool=false,
                 ifdivEtarg::Bool=false,ifcurlEtarg::Bool=false,
                 nd::Integer=1,thresh::Float64=1e-16)

    zk = complex(zk)
    
    # default values

    vals = FMMVals()
    
    zero = ComplexF64(0)
    zero3 = zeros(ComplexF64,3)
    
    Etarg = zero3
    divEtarg = zero
    curlEtarg = zero3

    ier = Integer(0)

    # check inputs

    ifE = false; ifdivE = false; ifcurlE = false

    anyfail, n, nt, ifA, ifB, iflam = (
        emfmm3dinputcheck(sources,A,B,lambda,
                          targets,ifE,ifdivE,ifcurlE,
                          ifEtarg,ifdivEtarg,
                          ifcurlEtarg,nd) )
    if anyfail
        return vals
    end

    if (ifA == 0); A = zero3 end
    if (ifB == 0); B = zero3 end
    if (iflam == 0); lambda = zero end    
    if (nt == 0); targets = zero3 end

    if (nt == 0)
        return vals
    end
    
    # allocate memory for return values

    ifE1 = Integer(0); ifcurlE1 = Integer(0)
    ifdivE1 = Integer(0)
    
    if ifEtarg
        if nd > 1
            Etarg = zeros(ComplexF64,nd,3,nt)
        else
            Etarg = zeros(ComplexF64,3,nt)
        end
        ifE1 = 1
    end

    if ifcurlEtarg
        if nd > 1
            curlEtarg = zeros(ComplexF64,nd,3,nt)
        else
            curlEtarg = zeros(ComplexF64,3,nt)
        end
        ifcurlE1 = 1
    end

    if ifdivEtarg
        if nd > 1
            divEtarg = zeros(ComplexF64,nd,nt)
        else
            divEtarg = zeros(ComplexF64,nt)
        end
        ifdivE1 = 1
    end

# actually call routine    
#
# subroutine em3ddirect(nd,zk,ns,source,ifa_vect,a_vect,&
#  ifb_vect,b_vect,iflambda,lambda,nt,targets,ifE,E,ifcurlE,curlE,&
#  ifdivE,divE,thresh)

    ccall((:em3ddirect_,libfmm3d),Cvoid,(Fi,Fc,Fi,Fd,Fi,Fc,
                                         Fi,Fc,Fi,Fc,Fi,Fd,
                                         Fi,Fc,Fi,Fc,Fi,Fc,
                                         Fd),
          nd,zk,n,sources,ifA,A,ifB,B,iflam,lambda,nt,targets,
          ifE1,Etarg,ifcurlE1,curlEtarg,ifdivE1,divEtarg,thresh)

    if ifEtarg; vals.Etarg = Etarg end
    if ifdivEtarg; vals.divEtarg = divEtarg end
    if ifcurlEtarg; vals.curlEtarg = curlEtarg end
    
    return vals

end



end # module
