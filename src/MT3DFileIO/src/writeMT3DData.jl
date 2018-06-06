"""
`writeMT3DData` writes data for 3D MT problem.

Input:
    datafile :: String    - the name of the data file.
    datInfo  :: DataInfo  - accompanying informations about data.
    predData :: Vector    - data values.
    dataErr  :: Vector    - data standard errors.

"""
function writeMT3DData(datafile::String, datInfo::DataInfo, predData::Array,
                       dataErr::Vector{Float64}=zeros(0))

    datID = open(datafile, "w")


    @printf(datID,"%-18s %s\n","# Format:","MT3DData_1.0")
    descrb  = "Data file generated at " * Libc.strftime(time())
    @printf(datID,"%-18s %s\n","# Description:", descrb)

    # receiver locations
    rxLoc = datInfo.rxLoc
    nr = size(rxLoc, 1)
    @printf(datID,"%-25s %4d\n","Receiver Location (m):", nr)
    @printf(datID,"%-10s %-12s %-12s %s\n","#","X","Y","Z")
    for i = 1:nr
        @printf(datID,"%12.2f %12.2f %12.2f\n",rxLoc[i,1], rxLoc[i,2], rxLoc[i,3])
    end

    # frequencies
    freqs = datInfo.freqs
    nF = length(freqs)
    @printf(datID,"%-20s %4d\n","Frequencies (Hz):",nF)
    for i = 1:nF
        @printf(datID,"%15.5e\n",freqs[i])
    end

    # data type, data components
    dataType = datInfo.dataType
    dataComp = datInfo.dataComp
    nDC = length(dataComp)
    @printf(datID,"DataType:  %s\n", dataType)
    @printf(datID,"DataComp:  %4d\n", nDC)
    for i = 1:nDC
        @printf(datID,"%s\n",dataComp[i])
    end


    if isempty(dataErr)
        dataErr = abs.(predData) * 0.05
    elseif length(dataErr) .== 1
        dataErr = abs.(predData) * dataErr[1]
    end

    rxID   = datInfo.rxID
    freqID = datInfo.freqID
    dcID   = datInfo.dcID
    nData  = size(predData,1)

    @printf(datID,"Data Block:  %6d\n", nData)

    if eltype(predData) <: Complex     # iseltype(predData, Complex)
        @printf(datID,"# %-8s %-7s %-9s %-15s %-15s %s\n","FreqNo.","RxNo.",
                "DCompNo.","RealValue","ImagValue","Error")
        for i = 1:nData
            @printf(datID,"%6d %8d %7d %15.6e %15.6e %15.6e\n", freqID[i], rxID[i],
                    dcID[i], real(predData[i]), imag(predData[i]), dataErr[i])
        end
    else
        @printf(datID,"# %-8s %-7s %-11s %-14s %s\n","FreqNo.","RxNo.",
                "DCompNo.","Value","Error")
        for i = 1:nData
            @printf(datID,"%6d %8d %7d %15.6e %15.6e\n", freqID[i], rxID[i],
                    dcID[i], predData[i], dataErr[i])
        end
    end

    close(datID)

end
