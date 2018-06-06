#
# TensorUtils defines some routines for Euler's rotation
# (c) HB, OUC, 28 April, 2016
#

#
# rotateTensor3!
#
#
function rotateTensor3!{T<:Real}(sigma::Array{T,2})

    # from degree to radian. strike, dip, slant
    radS = sigma[:,4] * pi/180
    radD = sigma[:,5] * pi/180
    radL = sigma[:,6] * pi/180

    sinS = sin.(radS)
    cosS = cos.(radS)
    sinD = sin.(radD)
    cosD = cos.(radD)
    sinL = sin.(radL)
    cosL = cos.(radL)

    ns = size(sigma,1)
    sigxx = zeros(Float64, ns)
    sigyy = zeros(Float64, ns)
    sigzz = zeros(Float64, ns)
    sigxy = zeros(Float64, ns)
    sigxz = zeros(Float64, ns)
    sigyz = zeros(Float64, ns)

    for j=1:ns
        rzS = [cosS[j]  -sinS[j]  0;
               sinS[j]   cosS[j]  0;
                 0         0      1]
        rxD = [1     0          0;
               0   cosD[j]  -sinD[j];
               0   sinD[j]   cosD[j]]
        rzL = [cosL[j]  -sinL[j]  0;
               sinL[j]   cosL[j]  0;
                 0         0      1]
        R = rzS * rxD * rzL
        sigM3 = spdiagm([sigma[j,1];sigma[j,2];sigma[j,3]])
        sigM3 = R * sigM3 * R'
        sigxx[j] = sigM3[1,1]
        sigxy[j] = sigM3[2,1]
        sigxz[j] = sigM3[3,1]
        sigyy[j] = sigM3[2,2]
        sigyz[j] = sigM3[3,2]
        sigzz[j] = sigM3[3,3]
    end

    sigma[:,1] = sigxx
    sigma[:,2] = sigyy
    sigma[:,3] = sigzz
    sigma[:,4] = sigxy
    sigma[:,5] = sigxz
    sigma[:,6] = sigyz

    return sigma

end  # rotateTensor3!



function rotateTensor3Inv!{T<:Real}(sigma::Array{T,2})
end  # rotateTensor3Inv!
