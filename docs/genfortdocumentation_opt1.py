from numpy import *

intro = "This subroutine evaluates the "
pgstr = ["potential ", "potential and it's gradient "]
intro2 = "\n\n  .. math::\n\n"
eq_start = "u(x) = \sum_{j=1}^{N} "
str1 = "c_{j} G_{k}(x-x_{j})"
str2 = "v_{j} \cdot \\nabla G_{k}(x-x_{j})"
str3 = str1+" + "+str2
eq_cjs = [eq_start+str1,eq_start+str2,eq_start+str3]

stflag = ["at the source locations $x=x_{j}$.", "at the target locations $x=t_{i}$.", "at the source and target locations $x=x_{j},t_{i}$."]
intro3 = "When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the sum."
inp_args ="Input arguments:"


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

c_opts = ['c','d','cd']
st_opts = ['stos','stot','stost']
p_opts = ['p','g']


f1 = open('fortrandocs_helm_opt1.raw','w')
for i in range(3):
    for j in range(3):
        for k in range(2):
            subname = 'subroutine hfmm3dpart'+st_opts[i]+c_opts[j]+p_opts[k]+'(eps,zk,nsource,source,'
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
            
            f1.writelines('.. code:: fortran\n\n')
            f1.writelines(sp1+subname+'\n\n')
            f1.writelines(intro+pgstr[k]+'\n'+intro2)
            f1.writelines("      "+eq_cjs[j]+"\n\n"+stflag[i]+" "+intro3+"\n\n")
            f1.writelines(inp_args+"\n\n")
            f1.writelines(sp1+eps_txt[0]+"\n"+sp2+eps_txt[1]+"\n")
            f1.writelines(sp1+zk_txt[0]+"\n"+sp2+zk_txt[1]+"\n")
            f1.writelines(sp1+ns_txt[0]+"\n"+sp2+ns_txt[1]+"\n")
            f1.writelines(sp1+src_txt[0]+"\n"+sp2+src_txt[1]+"\n")
            if(j==0):
                f1.writelines(sp1+charge_txt[0]+"\n"+sp2+charge_txt[1]+"\n")
            if(j==1):
                f1.writelines(sp1+dipole_txt[0]+"\n"+sp2+dipole_txt[1]+"\n")
            if(j==2):
                f1.writelines(sp1+charge_txt[0]+"\n"+sp2+charge_txt[1]+"\n")
                f1.writelines(sp1+dipole_txt[0]+"\n"+sp2+dipole_txt[1]+"\n")
            if(i>0):
                f1.writelines(sp1+nt_txt[0]+"\n"+sp2+nt_txt[1]+"\n")
                f1.writelines(sp1+targ_txt[0]+"\n"+sp2+targ_txt[1]+"\n")
            f1.writelines("\n\n")
            f1.writelines(inp_returns+"\n\n")
            if(i==0 or i==2):
                f1.writelines(sp1+pot_txt[0]+"\n"+sp2+pot_txt[1]+"\n")
                if(k==1):
                    f1.writelines(sp1+grad_txt[0]+"\n"+sp2+grad_txt[1]+"\n")
            if(i==1 or i==2):
                f1.writelines(sp1+pottarg_txt[0]+"\n"+sp2+pottarg_txt[1]+"\n")
                if(k==1):
                    f1.writelines(sp1+gradtarg_txt[0]+"\n"+sp2+gradtarg_txt[1]+"\n")
            f1.writelines('\n\n\n')

f1.close()

            
            
