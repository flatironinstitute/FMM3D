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

@testset "testing Stokes utilities" begin

    n = 1000
    nt = 1100
    sources = randn(3,n)
    targets = randn(3,nt)
    stoklet = randn(3,n)
    strslet = randn(3,n)
    strsvec = randn(3,n)

    vals1c = st3ddir(sources,targets,stoklet=stoklet,
                    ppregt=3,thresh=1e-16)

    vals1cs = st3ddir(sources,sources,stoklet=stoklet,
                    ppregt=3,thresh=1e-16)

    vals1d = st3ddir(sources,targets,strslet=strslet,strsvec=strsvec,
                    ppregt=3,thresh=1e-16)

    vals1ds = st3ddir(sources,sources,strslet=strslet,strsvec=strsvec,
                    ppregt=3,thresh=1e-16)

    vals1dc = st3ddir(sources,targets,stoklet=stoklet,
                     strslet=strslet,strsvec=strsvec,ppregt=3,thresh=1e-16)

    vals1dcs = st3ddir(sources,sources,stoklet=stoklet,
                     strslet=strslet,strsvec=strsvec,ppregt=3,thresh=1e-16)

    eps = 1e-12
    vals2c = stfmm3d(eps,sources,targets=targets,stoklet=stoklet,
                    ppreg=3,ppregt=3)
    vals2d = stfmm3d(eps,sources,targets=targets,strslet=strslet,strsvec=strsvec,
                    ppreg=3,ppregt=3)
    vals2dc = stfmm3d(eps,sources,targets=targets,stoklet=stoklet,
                     strslet=strslet,strsvec=strsvec,ppreg=3,ppregt=3)

    @test norm(vals2c.pot-vals1cs.pottarg)/norm(vals1cs.pottarg) < eps
    @test norm(vals2c.pre-vals1cs.pretarg)/norm(vals1cs.pretarg) < eps
    @test norm(vals2c.grad-vals1cs.gradtarg)/norm(vals1cs.gradtarg) < eps    
    @test norm(vals2c.pottarg-vals1c.pottarg)/norm(vals1c.pottarg) < eps
    @test norm(vals2c.gradtarg-vals1c.gradtarg)/norm(vals1c.gradtarg) < eps
    @test norm(vals2c.pretarg-vals1c.pretarg)/norm(vals1c.pretarg) < eps        

    @test norm(vals2d.pot-vals1ds.pottarg)/norm(vals1ds.pottarg) < eps
    @test norm(vals2d.pre-vals1ds.pretarg)/norm(vals1ds.pretarg) < eps
    @test norm(vals2d.grad-vals1ds.gradtarg)/norm(vals1ds.gradtarg) < eps
    @test norm(vals2d.pottarg-vals1d.pottarg)/norm(vals1d.pottarg) < eps
    @test norm(vals2d.pretarg-vals1d.pretarg)/norm(vals1d.pretarg) < eps
    @test norm(vals2d.gradtarg-vals1d.gradtarg)/norm(vals1d.gradtarg) < eps        

    @test norm(vals2dc.pot-vals1dcs.pottarg)/norm(vals1dcs.pottarg) < eps
    @test norm(vals2dc.pre-vals1dcs.pretarg)/norm(vals1dcs.pretarg) < eps
    @test norm(vals2dc.grad-vals1dcs.gradtarg)/norm(vals1dcs.gradtarg) < eps
    @test norm(vals2dc.pottarg-vals1dc.pottarg)/norm(vals1dc.pottarg) < eps
    @test norm(vals2dc.pretarg-vals1dc.pretarg)/norm(vals1dc.pretarg) < eps    
    @test norm(vals2dc.gradtarg-vals1dc.gradtarg)/norm(vals1dc.gradtarg) < eps    
    

end

@testset "testing electromagnetics utilities" begin

    n = 10
    nt = 11
    sources = randn(3,n)
    targets = randn(3,nt)
    A = randn(3,n) + im*randn(3,n)
    B = randn(3,n) + im*randn(3,n)
    lambda = randn(n) + im*randn(n)
    zk = 2.1

    vals1A = em3ddir(zk,sources,targets,A=A,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1As = em3ddir(zk,sources,sources,A=A,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1B = em3ddir(zk,sources,targets,B=B,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1Bs = em3ddir(zk,sources,sources,B=B,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1lambda = em3ddir(zk,sources,targets,lambda=lambda,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1lambdas = em3ddir(zk,sources,sources,lambda=lambda,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1all = em3ddir(zk,sources,targets,A=A,B=B,lambda=lambda,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    vals1alls = em3ddir(zk,sources,sources,A=A,B=B,lambda=lambda,
                     ifEtarg=true,ifdivEtarg=true,
                     ifcurlEtarg=true,thresh=1e-16)

    eps = 1e-12
    vals2A = emfmm3d(eps,zk,sources,targets=targets,A=A,
                    ifE=true,ifdivE=true,ifcurlE=true,
                    ifEtarg=true,ifdivEtarg=true,ifcurlEtarg=true)
    
    vals2B = emfmm3d(eps,zk,sources,targets=targets,B=B,
                    ifE=true,ifdivE=true,ifcurlE=true,
                    ifEtarg=true,ifdivEtarg=true,ifcurlEtarg=true)
    
    vals2lambda = emfmm3d(eps,zk,sources,targets=targets,lambda=lambda,
                    ifE=true,ifdivE=true,ifcurlE=true,
                    ifEtarg=true,ifdivEtarg=true,ifcurlEtarg=true)
    
    vals2all = emfmm3d(eps,zk,sources,targets=targets,A=A,B=B,lambda=lambda,
                    ifE=true,ifdivE=true,ifcurlE=true,
                    ifEtarg=true,ifdivEtarg=true,ifcurlEtarg=true)

    function absrelerr(vex,v)
        return norm(vex-v)/max(norm(vex),1)
    end
    
    @test absrelerr(vals1As.Etarg,vals2A.E) < eps
    @test absrelerr(vals1A.Etarg,vals2A.Etarg) < eps    
    @test absrelerr(vals1As.divEtarg,vals2A.divE) < eps
    @test absrelerr(vals1A.divEtarg,vals2A.divEtarg) < eps    
    @test absrelerr(vals1As.curlEtarg,vals2A.curlE) < eps
    @test absrelerr(vals1A.curlEtarg,vals2A.curlEtarg) < eps      
                  
    @test absrelerr(vals1Bs.Etarg,vals2B.E) < eps
    @test absrelerr(vals1B.Etarg,vals2B.Etarg) < eps    
    @test absrelerr(vals1Bs.divEtarg,vals2B.divE) < eps
    @test absrelerr(vals1B.divEtarg,vals2B.divEtarg) < eps    
    @test absrelerr(vals1Bs.curlEtarg,vals2B.curlE) < eps
    @test absrelerr(vals1B.curlEtarg,vals2B.curlEtarg) < eps      
                  
    @test absrelerr(vals1lambdas.Etarg,vals2lambda.E) < eps
    @test absrelerr(vals1lambda.Etarg,vals2lambda.Etarg) < eps    
    @test absrelerr(vals1lambdas.divEtarg,vals2lambda.divE) < eps
    @test absrelerr(vals1lambda.divEtarg,vals2lambda.divEtarg) < eps    
    @test absrelerr(vals1lambdas.curlEtarg,vals2lambda.curlE) < eps
    @test absrelerr(vals1lambda.curlEtarg,vals2lambda.curlEtarg) < eps      
                  
    @test absrelerr(vals1alls.Etarg,vals2all.E) < eps
    @test absrelerr(vals1all.Etarg,vals2all.Etarg) < eps    
    @test absrelerr(vals1alls.divEtarg,vals2all.divE) < eps
    @test absrelerr(vals1all.divEtarg,vals2all.divEtarg) < eps    
    @test absrelerr(vals1alls.curlEtarg,vals2all.curlE) < eps
    @test absrelerr(vals1all.curlEtarg,vals2all.curlEtarg) < eps      

end


@testset "testing lower level routines" begin


    # besseljs3d

    nterms = 10
    z = 1.1 + im*1.2
    scale = 1.3
    ifder = 1

    fj10 = (-3.5183264829466616750769540*(1e-9) +
        im*8.88237983492960538163695769*(1e-9))

    fjs, fjder = besseljs3d(nterms,z,scale=scale,ifder=ifder)

    @test (abs(fj10-fjs[11]*(scale^10))/abs(fj10) < 1e-10)

end
