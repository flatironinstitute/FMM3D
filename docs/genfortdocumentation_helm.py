from numpy import *

intro = "This subroutine evaluates the "
pgstr = ["potential ", "potential and its gradient "]
intro2 = "\n\n  .. math::\n\n"

eq_start = "u(x) = "
eq_start2 = "\sum_{j=1}^{N} "
str1 = "c_{j} \\frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}"
str2 = "v_{j} \cdot \\nabla \\left( \\frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\\right)" 
str3 = str1+" - "+str2
eq_cjs = [eq_start+eq_start2+str1,eq_start+"-"+eq_start2+str2,eq_start+eq_start2+str3]

eq_start = "u_{\ell}(x) = "
eq_start2 = "\sum_{j=1}^{N} "
str1 = "c_{\ell,j} \\frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}"
str2 = "v_{\ell,j} \cdot \\nabla \\left( \\frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\\right)" 
str3 = str1+" - "+str2
eq_nd_cjs = [eq_start+eq_start2+str1,eq_start+"-"+eq_start2+str2,eq_start+eq_start2+str3]

stflag = ["at the source locations $x=x_{j}$.", "at the target locations $x=t_{i}$.", "at the source and target locations $x=x_{j},t_{i}$."]
intro3 = "When $x=x_{j}$, the term corresponding to $x_{j}$ is "
intro3_cont = "dropped from the sum."
inp_args ="Input arguments:"


nd_txt = ["-    nd: integer","number of densities"]
eps_txt = ["-    eps: double precision","precision requested"]
zk_txt = ["-    zk: double complex","Helmholtz parameter, k"]
ns_txt = ["-    nsource: integer","Number of sources"]
src_txt =["-    source: double precision(3,nsource)","Source locations, $x_{j}$"]
charge_txt =["-    charge: double complex(nsource)","Charge strengths, $c_{j}$"]
charge_nd_txt =["-    charge: double complex(nd,nsource)","Charge strengths, $c_{\ell,j}$"]
dipole_txt =["-    dipvec: double complex(3,nsource)","Dipole strengths, $v_{j}$"]
dipole_nd_txt =["-    dipvec: double complex(nd,3,nsource)","Dipole strengths, $v_{\ell,j}$"]
nt_txt = ["-    ntarg: integer","Number of targets"]
targ_txt =["-    targ: double precision(3,ntarg)","Target locations, $t_{i}$"]

inp_returns = "Output arguments:"
pot_txt =["-    pot: double complex(nsource)","Potential at source locations, $u(x_{j})$"]
grad_txt =["-    grad: double complex(3,nsource)","Gradient at source locations, $\\nabla u(x_{j})$"]
pottarg_txt =["-    pottarg: double complex(ntarg)","Potential at target locations, $u(t_{i})$"]
gradtarg_txt =["-    gradtarg: double complex(3,ntarg)","Gradient at target locations, $\\nabla u(t_{i})$"]

pot_nd_txt =["-    pot: double complex(nd,nsource)","Potential at source locations, $u_{\ell}(x_{j})$"]
grad_nd_txt =["-    grad: double complex(nd,3,nsource)","Gradient at source locations, $\\nabla u_{\ell}(x_{j})$"]
pottarg_nd_txt =["-    pottarg: double complex(nd,ntarg)","Potential at target locations, $u_{\ell}(t_{i})$"]
gradtarg_nd_txt =["-    gradtarg: double complex(nd,3,ntarg)","Gradient at target locations, $\\nabla u_{\ell}(t_{i})$"]

ier_txt= ["-    ier: integer","Error flag; ier=0 implies successful execution, and ier=4/8 implies insufficient memory"] 
sp1 = "  "
sp1 = "  "
sp2 = "          "

c_opts = ['_c','_d','_cd']
c_opts2 = ['Charges','Dipoles', 'Charges and Dipoles']
st_opts = ['_s','_t','_st']
st_opts2 = ['Sources','Targets','Sources and Targets']
p_opts = ['_p','_g']
p_opts2 = ['Potential','Potential and Gradient']

###
#  Generate webpage for helmholtz
#
f1 = open('fortrandocs_helm.raw','w')
for i in range(3):
    for j in range(3):
        for k in range(2):
            subname0 = 'subroutine hfmm3d'+st_opts[i]+c_opts[j]+p_opts[k]+'(eps,zk,nsource,source,'
            subname1 = 'subroutine hfmm3d'+st_opts[i]+c_opts[j]+p_opts[k]+'_vec(nd,eps,zk,nsource,source,'
            subname = ''
            if(j==0 or j==2):
                subname=subname+'charge,'
            if(j==1 or j==2):
                subname=subname+'dipvec,'
            if(i == 2 and k==0):
                subname=subname+'pot,'
            if(i == 2 and k==1):
                subname=subname+'pot,grad,'
            if(i>0):
                subname=subname+'ntarg,targ,'
            if(i==0 and k==0):
                subname=subname+'pot,ier)'
            if(i==0 and k==1):
                subname=subname+'pot,grad,ier)'
            if(i>=1 and k==0):
                subname=subname+'pottarg,ier)'
            if(i>=1 and k==1):
                subname=subname+'pottarg,gradtarg,ier)'

            
            str_ini = 'h'+st_opts[i][1::]+c_opts[j][1::]+p_opts[k][1::]

            subname3 = 'hfmm3d'+st_opts[i]+c_opts[j]+p_opts[k]
            f1.writelines('.. _'+str_ini+':\n\n')
            f1.writelines(subname3+'\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n')
            f1.writelines('- Evaluation points: '+st_opts2[i]+'\n')
            f1.writelines('- Interaction kernel: '+c_opts2[j]+'\n')
            f1.writelines('- Outputs requested: '+p_opts2[k]+'\n\n')
            f1.writelines('-------------------------------------\n\n')
            f1.writelines('.. code:: fortran\n\n')
            f1.writelines(sp1+subname0+subname+'\n\n')
            f1.writelines(intro+pgstr[k]+'\n'+intro2)
            f1.writelines("      "+eq_cjs[j]+"\n\n"+stflag[i]+" "+intro3+intro3_cont+"\n\n")
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
            f1.writelines(sp1+ier_txt[0]+"\n"+sp2+ier_txt[1]+"  \n")
            f1.writelines("\n\n--------------------------------\n\nVectorized version: \n\n")
            f1.writelines('.. code:: fortran\n\n')
            f1.writelines(sp1+subname1+subname+'\n\n')
            f1.writelines(intro+pgstr[k]+'\n'+intro2)
            f1.writelines("      "+eq_nd_cjs[j]+"\n\n"+stflag[i]+" "+intro3+intro3_cont+"\n\n")
            f1.writelines(inp_args+"\n\n")
            f1.writelines(sp1+nd_txt[0]+"\n"+sp2+nd_txt[1]+"\n")
            if(j==0):
                f1.writelines(sp1+charge_nd_txt[0]+"\n"+sp2+charge_nd_txt[1]+"\n")
            if(j==1):
                f1.writelines(sp1+dipole_nd_txt[0]+"\n"+sp2+dipole_nd_txt[1]+"\n")
            if(j==2):
                f1.writelines(sp1+charge_nd_txt[0]+"\n"+sp2+charge_nd_txt[1]+"\n")
                f1.writelines(sp1+dipole_nd_txt[0]+"\n"+sp2+dipole_nd_txt[1]+"\n")
            f1.writelines("\n\n")
            f1.writelines(inp_returns+"\n\n")
            if(i==0 or i==2):
                f1.writelines(sp1+pot_nd_txt[0]+"\n"+sp2+pot_nd_txt[1]+"\n")
                if(k==1):
                    f1.writelines(sp1+grad_nd_txt[0]+"\n"+sp2+grad_nd_txt[1]+"\n")
            if(i==1 or i==2):
                f1.writelines(sp1+pottarg_nd_txt[0]+"\n"+sp2+pottarg_nd_txt[1]+"\n")
                if(k==1):
                    f1.writelines(sp1+gradtarg_nd_txt[0]+"\n"+sp2+gradtarg_nd_txt[1]+"\n")
            f1.writelines(sp1+ier_txt[0]+"\n"+sp2+ier_txt[1]+"  \n")
            f1.writelines("\n\n.. container:: rttext\n\n  `Back to Helmholtz FMM <fortran-c.html#helm>`__")
            f1.writelines("\n\n.. container:: rttext\n\n  `Back to top <fortran-c.html#fcexmp>`__\n\n\n")

f1.close()

######
#  Write fortran comments for helmholtz
#
eq_start = "u(x) = "
eq_start2 = "\sum_{j=1}^{N} "
str1 = "c_{j} \\frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}"
str2 = "v_{j} \cdot \\nabla \\left( \nc            \\frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\\right)" 
str3 = str1+" - \nc            "+str2
eq_cjs_fort = [eq_start+eq_start2+str1,eq_start+"-"+eq_start2+str2,eq_start+eq_start2+str3]

eq_start = "u_{\ell}(x) = "
eq_start2 = "\sum_{j=1}^{N} "
str1 = "c_{\ell,j}\nc        \\frac{e^{ik\|x- x_{j}\|}}{\|x-x_{j}\|}"
str2 = "v_{\ell,j} \cdot \\nabla \\left( \nc        \\frac{e^{ik\|x-x_{j}\|}}{\|x-x_{j}\|}\\right)" 
str3 = str1+" - \nc            "+str2
eq_cjs_fort_nd = [eq_start+eq_start2+str1,eq_start+"-"+eq_start2+str2,eq_start+eq_start2+str3]


f1 = open('fortrandocs_helm2.raw','w')
f2 = open('fortrandocs_helm2_vec.raw','w')
f3 = open('fortrandocs_helm_header.raw','w')
f4 = open('fortrandocs_helm_header_vec.raw','w')
for i in range(3):
    for j in range(3):
        for k in range(2):
            f1.writelines('c-------------------------------------\n')

            f3.writelines('c  -hfmm3d'+st_opts[i]+c_opts[j]+p_opts[k]+'\n')
            f3.writelines('c    - Evaluation points: '+st_opts2[i]+'\n')
            f3.writelines('c    - Interaction kernel: '+c_opts2[j]+'\n')
            f3.writelines('c    - Outputs requested: '+p_opts2[k]+'\n')
            f3.writelines('c-------------------------------------\n')

            f4.writelines('c  -hfmm3d'+st_opts[i]+c_opts[j]+p_opts[k]+'_vec\n')
            f4.writelines('c    - Evaluation points: '+st_opts2[i]+'\n')
            f4.writelines('c    - Interaction kernel: '+c_opts2[j]+'\n')
            f4.writelines('c    - Outputs requested: '+p_opts2[k]+'\n')
            f4.writelines('c-------------------------------------\n')

            f1.writelines('c\n')
            f1.writelines("c  "+intro+pgstr[k]+'\n')
            f1.writelines("c      "+eq_cjs_fort[j]+"\nc\nc  "+stflag[i]+"\nc  "+intro3+"\nc  "+intro3_cont+"\nc\n")
            f1.writelines("c  "+inp_args+"\nc\n")
            f1.writelines("c  "+sp1+eps_txt[0]+"\nc"+sp2+eps_txt[1]+"\n")
            f1.writelines("c  "+sp1+zk_txt[0]+"\nc"+sp2+zk_txt[1]+"\n")
            f1.writelines("c  "+sp1+ns_txt[0]+"\nc"+sp2+ns_txt[1]+"\n")
            f1.writelines("c  "+sp1+src_txt[0]+"\nc"+sp2+src_txt[1]+"\n")
            if(j==0):
                f1.writelines("c  "+sp1+charge_txt[0]+"\nc"+sp2+charge_txt[1]+"\n")
            if(j==1):
                f1.writelines("c  "+sp1+dipole_txt[0]+"\nc"+sp2+dipole_txt[1]+"\n")
            if(j==2):
                f1.writelines("c  "+sp1+charge_txt[0]+"\nc"+sp2+charge_txt[1]+"\n")
                f1.writelines("c  "+sp1+dipole_txt[0]+"\nc"+sp2+dipole_txt[1]+"\n")
            if(i>0):
                f1.writelines("c  "+sp1+nt_txt[0]+"\nc"+sp2+nt_txt[1]+"\n")
                f1.writelines("c  "+sp1+targ_txt[0]+"\nc"+sp2+targ_txt[1]+"\n")
            f1.writelines("c\nc\n")
            f1.writelines("c  "+inp_returns+"\nc\n")
            if(i==0 or i==2):
                f1.writelines("c  "+sp1+pot_txt[0]+"\nc"+sp2+pot_txt[1]+"\nc")
                if(k==1):
                    f1.writelines("  "+sp1+grad_txt[0]+"\nc"+sp2+grad_txt[1]+"\nc")
            if(i==1 or i==2):
                f1.writelines("c  "+sp1+pottarg_txt[0]+"\nc"+sp2+pottarg_txt[1]+"\nc")
                if(k==1):
                    f1.writelines("  "+sp1+gradtarg_txt[0]+"\nc"+sp2+gradtarg_txt[1]+"\nc")
            f1.writelines("  "+sp1+ier_txt[0]+"\nc"+sp2+ier_txt[1]+"  \nc")
            f1.writelines("\nc\nc--------------------------------\nc\n")

            f2.writelines('c-------------------------------------\n')
            f2.writelines('c\n')
            f2.writelines("c  "+intro+pgstr[k]+'\n')
            f2.writelines("c      "+eq_cjs_fort_nd[j]+"\nc\nc  "+stflag[i]+"\nc  "+intro3+"\nc  "+intro3_cont+"\nc\n")
            f2.writelines("c  "+inp_args+"\nc\n")
            f2.writelines("c  "+sp1+nd_txt[0]+"\nc"+sp2+nd_txt[1]+"\n")
            f2.writelines("c  "+sp1+eps_txt[0]+"\nc"+sp2+eps_txt[1]+"\n")
            f2.writelines("c  "+sp1+zk_txt[0]+"\nc"+sp2+zk_txt[1]+"\n")
            f2.writelines("c  "+sp1+ns_txt[0]+"\nc"+sp2+ns_txt[1]+"\n")
            f2.writelines("c  "+sp1+src_txt[0]+"\nc"+sp2+src_txt[1]+"\n")
            if(j==0):
                f2.writelines("c  "+sp1+charge_nd_txt[0]+"\nc"+sp2+charge_nd_txt[1]+"\n")
            if(j==1):
                f2.writelines("c  "+sp1+dipole_nd_txt[0]+"\nc"+sp2+dipole_nd_txt[1]+"\n")
            if(j==2):
                f2.writelines("c  "+sp1+charge_nd_txt[0]+"\nc"+sp2+charge_nd_txt[1]+"\n")
                f2.writelines("c  "+sp1+dipole_nd_txt[0]+"\nc"+sp2+dipole_nd_txt[1]+"\n")
            if(i>0):
                f2.writelines("c  "+sp1+nt_txt[0]+"\nc"+sp2+nt_txt[1]+"\n")
                f2.writelines("c  "+sp1+targ_txt[0]+"\nc"+sp2+targ_txt[1]+"\n")
            f2.writelines("c\nc\n")
            f2.writelines("c  "+inp_returns+"\nc\n")
            if(i==0 or i==2):
                f2.writelines("c  "+sp1+pot_nd_txt[0]+"\nc"+sp2+pot_nd_txt[1]+"\nc")
                if(k==1):
                    f2.writelines("  "+sp1+grad_nd_txt[0]+"\nc"+sp2+grad_nd_txt[1]+"\nc")
            if(i==1 or i==2):
                f2.writelines("c  "+sp1+pottarg_nd_txt[0]+"\nc"+sp2+pottarg_nd_txt[1]+"\nc")
                if(k==1):
                    f2.writelines("  "+sp1+gradtarg_nd_txt[0]+"\nc"+sp2+gradtarg_nd_txt[1]+"\nc")
            f2.writelines("  "+sp1+ier_txt[0]+"\n"+sp2+ier_txt[1]+"  \nc")
            f2.writelines("\nc\nc--------------------------------\nc\n")
f1.close()
f2.close()
f3.close()
f4.close()
