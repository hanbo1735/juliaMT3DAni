"""
`writeMT3DResp` writes forward responses for 3D MT problem.

Input:
    respfile :: String    - the name of the response file.
    datInfo  :: DataInfo  - accompanying informations about data.
    fwdResp  :: Array     - forward responses.

"""
function writeMT3DResp{T<:Union{Float64, Complex128}}(respfile::String, datInfo::DataInfo,
                       fwdResp::Array{T, 3})

    fid = open(respfile, "w")

    @printf(fid,"%-18s %s\n","# Format:","MT3DResp_1.0")
    descrb  = "Response file generated at " * Libc.strftime(time())
    @printf(fid,"%-18s %s\n","# Description:", descrb)

    # receiver location
    rxLoc = datInfo.rxLoc
    nRx = size(rxLoc, 1)
    @printf(fid,"%-25s %4d\n","Receiver Location (m):", nRx)
    @printf(fid,"%-10s %-12s %-12s %s\n","#","X","Y","Z")
    for i = 1:nRx
        @printf(fid,"%12.2f %12.2f %12.2f\n", rxLoc[i,1], rxLoc[i,2], rxLoc[i,3])
    end

    # frequencies
    freqs = datInfo.freqs
    nFreq = length(freqs)
    @printf(fid,"%-20s %4d\n","Frequencies (Hz):", nFreq)
    for i = 1:nFreq
        @printf(fid,"%15.5e\n", freqs[i])
    end

    # data type, data components
    dataType = datInfo.dataType
    dataComp = datInfo.dataComp
    nDC = length(dataComp)
    @printf(fid,"DataType:  %s\n", dataType)
    @printf(fid,"DataComp:  %4d\n", nDC)
    for i = 1:nDC
        @printf(fid,"%s\n",dataComp[i])
    end


    nData  = nFreq * nRx

    @printf(fid,"%-15s %d\n","Data Block:",nData)


    if dataType == "Impedance"
        @printf(fid,"# %-8s %-9s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %s\n","FreqNo.","RxNo.",
                "RealZXX","ImagZXX", "RealZXY","ImagZXY","RealZYX","ImagZYX","RealZYY","ImagZYY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", j, k,
                    real(fwdResp[k,1,j]), imag(fwdResp[k,1,j]), real(fwdResp[k,2,j]), imag(fwdResp[k,2,j]),
                    real(fwdResp[k,3,j]), imag(fwdResp[k,3,j]), real(fwdResp[k,4,j]), imag(fwdResp[k,4,j]))
        end

    elseif dataType == "Impedance_Tipper"
        @printf(fid,"# %-8s %-9s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %s\n",
                "FreqNo.","RxNo.","RealZXX","ImagZXX", "RealZXY","ImagZXY","RealZYX","ImagZYX",
                "RealZYY","ImagZYY","RealTZX","ImagTZX","RealTZY","ImagTZY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n", j, k,
                    real(fwdResp[k,1,j]), imag(fwdResp[k,1,j]), real(fwdResp[k,2,j]), imag(fwdResp[k,2,j]),
                    real(fwdResp[k,3,j]), imag(fwdResp[k,3,j]), real(fwdResp[k,4,j]), imag(fwdResp[k,4,j]),
                    real(fwdResp[k,5,j]), imag(fwdResp[k,5,j]), real(fwdResp[k,6,j]), imag(fwdResp[k,6,j]))
        end

    elseif dataType == "Rho_Phs"
        @printf(fid,"# %-8s %-12s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %s\n","FreqNo.","RxNo.",
                "RhoXX","PhsXX","RhoXY","PhsXY","RhoYX","PhsYX","RhoYY","PhsYY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f\n", j, k,
                    fwdResp[k,1,j], fwdResp[k,2,j], fwdResp[k,3,j], fwdResp[k,4,j],
                    fwdResp[k,5,j], fwdResp[k,6,j], fwdResp[k,7,j], fwdResp[k,8,j])
        end

    elseif dataType == "Rho_Phs_Tipper"
        @printf(fid,"# %-8s %-12s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-13s %-15s %-15s %-15s %s\n","FreqNo.","RxNo.",
                "RhoXX","PhsXX","RhoXY","PhsXY","RhoYX","PhsYX","RhoYY","PhsYY",
                "RealTZX","ImagTZX","RealTZY","ImagTZY")

        for j in 1:nFreq, k in 1:nRx
            @printf(fid,"%6d %8d %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %13.2f %15.6e %15.6e %15.6e %15.6e\n", j, k,
                    fwdResp[k,1,j], fwdResp[k,2,j], fwdResp[k,3,j], fwdResp[k,4,j],
                    fwdResp[k,5,j], fwdResp[k,6,j], fwdResp[k,7,j], fwdResp[k,8,j],
                    fwdResp[k,9,j], fwdResp[k,10,j], fwdResp[k,11,j], fwdResp[k,12,j])
        end

    end  # dataType==?

    close(fid)

end
