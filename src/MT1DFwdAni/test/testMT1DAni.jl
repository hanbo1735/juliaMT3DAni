#
# script: test function in module MTFwdSolver
#
workspace()
using MT3DFwdSolver
import YeeGridOperators: CondTensor, EMTensorMesh
import Utilities.rotateTensor3!


#
println("=== test MT 1D anisotropic forward functions")
#

# The 1D Anisotropic model used in Pek's paper (CG, 2002)
sigma = [10000.  10000.  10000.    0.  0.  0.;
           200.  20000.    200.   15.  0.  0.;
          1000.   2000.   1000.  -75.  0.  0.;
           100.    100.    100.    0.  0.  0.]

sigma[:,1:3] = 1./sigma[:,1:3]
rotateTensor3!(sigma)
h = [0; 10.0; 18; 100]
zNode = cumsum(h) * 1000
freqs = [ 0.0001; 0.001; 0.01; 0.1; 1.; 10.; 100.]

#=
println("testing mt1DAniAnalyticField ...")
zNodeP = collect([zNode; zNode[end]+100])  # make sure length(zNode) > size(sigma,1)
nFreq = length(freqs)
nLayer = length(zNode)

Ex = zeros(Complex128, nLayer, nFreq)
Ey = copy(Ex)
Ez = zeros(Complex128, nLayer-1, nFreq)
Hx = copy(Ex)
Hy = copy(Ex)

eTop = [1.0+0im; 0]

for j=1:nFreq
    Ex[:,j], Ey[:,j], Ez[:,j], Hx[:,j], Hy[:,j] = mt1DAniAnalyticField(freqs[j], sigma, zNodeP, eTop; compH=true)
end


fid = open("field.txt", "w")

for iFreq=1:nFreq
    ex = Ex[:, iFreq]
    ey = Ey[:, iFreq]
    ez = [Ez[:, iFreq]; 0]
    hx = Hx[:, iFreq]
    hy = Hy[:, iFreq]

    for j=1:nLayer
    @printf(fid,"%4d %10.3f %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n",
            iFreq, zNode[j], real(ex[j]), imag(ex[j]), real(ey[j]), imag(ey[j]),
            real(ez[j]), imag(ez[j]), real(hx[j]), imag(hx[j]), real(hy[j]), imag(hy[j]))
    end
    @printf(fid,"\n")
end

close(fid)
=#

println("testing mt1DAniFwdPek ...")
Z1, Rho1, Phz1 = mt1DAniFwdPek(freqs, sigma, zNode)


println("testing mt1DAniFwd ...")
Z2, Rho2, Phz2 = mt1DAniFwd(freqs, sigma, zNode)

# Cross-check
println("|Z1-Z2|/|Z1|  |Rho1-Rho2|/|Rho1|  |Phz1-Phz2|/|Phz1|")
@printf("%6.2e %6.2e %6.2e\n",norm(Z1-Z2)/norm(Z1), norm(Rho1-Rho2)/norm(Rho1), norm(Phz1-Phz2)/norm(Phz1) )
println("")

println("=== All test passed.")
