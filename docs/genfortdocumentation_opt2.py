from numpy import *

sp2 = "          "

c_opts = ['c','d','cd']
st_opts = ['stos','stot','stost']
p_opts = ['p','g']

sttext = 'Evaluation points: '
stflag = ['sources; ','targets; ',' sources and targets; ']

cdtext = 'Interaction: '
cdflag = ['charges; ','dipoles; ','charges and dipoles; ']

pgtext = 'Output: '
pgflag = ['potential.','potential and gradient.']

sp1 = '    '

f1 = open('fortrandocs_helm_opt2.raw','w')
f1.write('List of interfaces\n\n')
for i in range(3):
    for j in range(3):
        for k in range(2):
            f1.writelines(sp1+'-  '+sttext+stflag[i]+cdtext+cdflag[j]+pgtext+pgflag[k]+'\n\n')
            subname0 = 'subroutine hfmm3dpart'+st_opts[i]+c_opts[j]+p_opts[k]+'(eps,zk,nsource,source,'
            subname1 = 'subroutine hfmm3dpart'+st_opts[i]+c_opts[j]+p_opts[k]+'_vec(nd,eps,zk,nsource,source,'
            subname=''
            if(j==0 or j==2):
                subname=subname+'charge,'
            if(j==1 or j==2):
                subname=subname+'dipvec,'
            if(i>0):
                subname=subname+'ntarg,targ,'
            if(i==0 and k==0):
                subname=subname+'pot)'
            if(i==0 and k==1):
                subname=subname+'pot,grad)'
            if(i==1 and k==0):
                subname=subname+'pottarg)'
            if(i==1 and k==1):
                subname=subname+'pottarg,gradtarg)'
            if(i==2 and k==0):
                subname=subname+'pot,pottarg)'
            if(i==2 and k==1):
                subname=subname+'pot,grad,pottarg,gradtarg)'
             
            f1.writelines(sp1+'.. code:: fortran\n\n')
            f1.writelines(sp1+sp1+subname0+subname+'\n\n')
            f1.writelines(sp1+'For vectorized routines use: \n\n')
            f1.writelines(sp1+'.. code:: fortran\n\n')
            f1.writelines(sp1+sp1+subname1+subname+'\n\n')


            
inp_args ="Input arguments:"


nd_txt = ["-    nd: integer","number of charge/dipole densities"]
eps_txt = ["-    eps: double precision","precision requested"]
zk_txt = ["-    zk: double complex","Helmholtz parameter (k)"]
ns_txt = ["-    nsource: integer","Number of sources (nsource)"]
src_txt =["-    source: double precision(3,nsource)","Source locations ($x_{j}$)"]
charge_txt =["-    charge: double complex(nsource)","Charge strengths ($c_{j}$)"]
dipole_txt =["-    dipvec: double complex(3,nsource)","Dipole strengths ($v_{j}$)"]
nt_txt = ["-    ntarg: integer","Number of targets"]
targ_txt =["-    targ: double precision(3,ntarg)","Target locations ($t_{i}$)"]

inp_returns = "Output arguments:"
pot_txt =["-    pot: double complex(nsource)","Potential at source locations ($u(x_{j})$)"]
grad_txt =["-    grad: double complex(3,nsource)","Gradient at source locations ($\\nabla u(x_{j})$)"]
pottarg_txt =["-    pottarg: double complex(ntarg)","Potential at target locations ($u(t_{i})$)"]
gradtarg_txt =["-    gradtarg: double complex(3,ntarg)","Gradient at target locations ($\\nabla u(t_{i})$)"]


sp1 = "  "
sp2 = "          "

f1.writelines(inp_args+"\n************************\n\n")
f1.writelines(sp1+nd_txt[0]+"\n"+sp2+nd_txt[1]+"\n")
f1.writelines(sp1+eps_txt[0]+"\n"+sp2+eps_txt[1]+"\n")
f1.writelines(sp1+zk_txt[0]+"\n"+sp2+zk_txt[1]+"\n")
f1.writelines(sp1+ns_txt[0]+"\n"+sp2+ns_txt[1]+"\n")
f1.writelines(sp1+src_txt[0]+"\n"+sp2+src_txt[1]+"\n")
f1.writelines(sp1+charge_txt[0]+"\n"+sp2+charge_txt[1]+"\n")
f1.writelines(sp1+dipole_txt[0]+"\n"+sp2+dipole_txt[1]+"\n")
f1.writelines(sp1+nt_txt[0]+"\n"+sp2+nt_txt[1]+"\n")
f1.writelines(sp1+targ_txt[0]+"\n"+sp2+targ_txt[1]+"\n")

f1.writelines("\n\n")
f1.writelines(inp_returns+"\n***********************\n")
f1.writelines(sp1+pot_txt[0]+"\n"+sp2+pot_txt[1]+"\n")
f1.writelines(sp1+grad_txt[0]+"\n"+sp2+grad_txt[1]+"\n")
f1.writelines(sp1+pottarg_txt[0]+"\n"+sp2+pottarg_txt[1]+"\n")
f1.writelines(sp1+gradtarg_txt[0]+"\n"+sp2+gradtarg_txt[1]+"\n")
f1.close()
