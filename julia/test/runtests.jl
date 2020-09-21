using FMM3D
using Test
using Random
using LinearAlgebra

@testset "testing Helmholtz utilities" begin

    n = 1000
    nt = 1100
    sources = randn(3,n)
    targets = randn(3,nt)
    charges = randn(n) + im*randn(n)
    dipvecs = randn(3,n) + im*randn(3,n)
    zk = 2.1
    eps = 1e-9

    vals1c = h3ddir(zk,sources,targets,charges=charges,
                    pgt=2,thresh=1e-16)

    vals1cs = h3ddir(zk,sources,sources,charges=charges,
                    pgt=2,thresh=1e-16)

    vals1d = h3ddir(zk,sources,targets,dipvecs=dipvecs,
                    pgt=2,thresh=1e-16)

    vals1ds = h3ddir(zk,sources,sources,dipvecs=dipvecs,
                    pgt=2,thresh=1e-16)

    vals1dc = h3ddir(zk,sources,targets,charges=charges,
                     dipvecs=dipvecs,pgt=2,thresh=1e-16)

    vals1dcs = h3ddir(zk,sources,sources,charges=charges,
                     dipvecs=dipvecs,pgt=2,thresh=1e-16)

    eps = 1e-12
    vals2c = hfmm3d(eps,zk,sources,targets=targets,charges=charges,
                    pg=2,pgt=2)
    vals2d = hfmm3d(eps,zk,sources,targets=targets,dipvecs=dipvecs,
                    pg=2,pgt=2)
    vals2dc = hfmm3d(eps,zk,sources,targets=targets,charges=charges,
                     dipvecs=dipvecs,pg=2,pgt=2)

    @test norm(vals2c.pot-vals1cs.pottarg)/norm(vals1cs.pottarg) < eps
    @test norm(vals2c.grad-vals1cs.gradtarg)/norm(vals1cs.gradtarg) < eps
    @test norm(vals2c.pottarg-vals1c.pottarg)/norm(vals1c.pottarg) < eps
    @test norm(vals2c.gradtarg-vals1c.gradtarg)/norm(vals1c.gradtarg) < eps    

    @test norm(vals2d.pot-vals1ds.pottarg)/norm(vals1ds.pottarg) < eps
    @test norm(vals2d.grad-vals1ds.gradtarg)/norm(vals1ds.gradtarg) < eps
    @test norm(vals2d.pottarg-vals1d.pottarg)/norm(vals1d.pottarg) < eps
    @test norm(vals2d.gradtarg-vals1d.gradtarg)/norm(vals1d.gradtarg) < eps    

    @test norm(vals2dc.pot-vals1dcs.pottarg)/norm(vals1dcs.pottarg) < eps
    @test norm(vals2dc.grad-vals1dcs.gradtarg)/norm(vals1dcs.gradtarg) < eps
    @test norm(vals2dc.pottarg-vals1dc.pottarg)/norm(vals1dc.pottarg) < eps
    @test norm(vals2dc.gradtarg-vals1dc.gradtarg)/norm(vals1dc.gradtarg) < eps    
    

end

@testset "testing Laplace utilities" begin

    n = 1000
    nt = 1100
    sources = randn(3,n)
    targets = randn(3,nt)
    charges = randn(n)
    dipvecs = randn(3,n)
    eps = 1e-9

    vals1c = l3ddir(sources,targets,charges=charges,
                    pgt=2,thresh=1e-16)

    vals1cs = l3ddir(sources,sources,charges=charges,
                    pgt=2,thresh=1e-16)

    vals1d = l3ddir(sources,targets,dipvecs=dipvecs,
                    pgt=2,thresh=1e-16)

    vals1ds = l3ddir(sources,sources,dipvecs=dipvecs,
                    pgt=2,thresh=1e-16)

    vals1dc = l3ddir(sources,targets,charges=charges,
                     dipvecs=dipvecs,pgt=2,thresh=1e-16)

    vals1dcs = l3ddir(sources,sources,charges=charges,
                     dipvecs=dipvecs,pgt=2,thresh=1e-16)

    eps = 1e-12
    vals2c = lfmm3d(eps,sources,targets=targets,charges=charges,
                    pg=2,pgt=2)
    vals2d = lfmm3d(eps,sources,targets=targets,dipvecs=dipvecs,
                    pg=2,pgt=2)
    vals2dc = lfmm3d(eps,sources,targets=targets,charges=charges,
                     dipvecs=dipvecs,pg=2,pgt=2)

    @test norm(vals2c.pot-vals1cs.pottarg)/norm(vals1cs.pottarg) < eps
    @test norm(vals2c.grad-vals1cs.gradtarg)/norm(vals1cs.gradtarg) < eps
    @test norm(vals2c.pottarg-vals1c.pottarg)/norm(vals1c.pottarg) < eps
    @test norm(vals2c.gradtarg-vals1c.gradtarg)/norm(vals1c.gradtarg) < eps    

    @test norm(vals2d.pot-vals1ds.pottarg)/norm(vals1ds.pottarg) < eps
    @test norm(vals2d.grad-vals1ds.gradtarg)/norm(vals1ds.gradtarg) < eps
    @test norm(vals2d.pottarg-vals1d.pottarg)/norm(vals1d.pottarg) < eps
    @test norm(vals2d.gradtarg-vals1d.gradtarg)/norm(vals1d.gradtarg) < eps    

    @test norm(vals2dc.pot-vals1dcs.pottarg)/norm(vals1dcs.pottarg) < eps
    @test norm(vals2dc.grad-vals1dcs.gradtarg)/norm(vals1dcs.gradtarg) < eps
    @test norm(vals2dc.pottarg-vals1dc.pottarg)/norm(vals1dc.pottarg) < eps
    @test norm(vals2dc.gradtarg-vals1dc.gradtarg)/norm(vals1dc.gradtarg) < eps    
    

end

