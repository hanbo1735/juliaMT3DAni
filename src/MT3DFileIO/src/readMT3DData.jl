"""
`readMT3DData` reads the data file for 3D MT problem.

Input:
    datafile::String   - the name of the data file.

Output:
    datInfo::DataInfo  - accompanying informations about data.
    obsData::Vector    - data values.
    dataErr::Vector    - data standard errors.

"""
function readMT3DData(datafile::String)

    if isfile(datafile)
        fid = open(datafile, "r")
    else
        error("$(datafile) does not exist, please try again!")
    end

    # declare variables to be output
    rxLoc = []
    freqs = []
    dataType = []
    dataComp = []
    rxID   = []
    freqID = []
    dcID   = []

    obsData = []
    dataErr = []


    while !eof(fid)

        cline = strip(readline(fid))

        # ignore empty lines and comments (preceded with #)
        isempty(cline)  && continue
        cline[1] == '#' && continue

        # data format, won't be used in this version
        if contains(cline, "Format")
            tmp = split(cline)
            format = tmp[2]

        # receiver locations
        elseif contains(cline, "Receiver Location")
            tmp = split(cline)
            nr  = parse(Int, tmp[end])
            nd  = 0
            rxLoc = zeros(nr, 3)
            while nd < nr
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                for j = 1:3
                    rxLoc[nd,j] = float(cline[j])
                end
            end

        # frequencies
        elseif contains(cline, "Frequencies")
            tmp = split(cline)
            nf  = parse(Int, tmp[end])
            nd  = 0
            freqs = zeros(nf)
            while nd < nf
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                #cline = split(cline)
                nd = nd + 1
                freqs[nd] = float(cline)
            end

        # data type, allowed values: Impedance, Impedance_Tipper,
        #                            Rho_Phs, Rho_Phs_Tipper
        elseif contains(cline, "DataType")
            tmp = split(cline)
            dataType = tmp[end]
            if !contains(dataType, "Impedance") && !contains(dataType, "Rho_Phs")
                error("Data type: $(dataType) is not supported!")
            end

            if contains(dataType, "Impedance")
                isComplex = true
            else
                isComplex = false
            end

        # data components, allowed values for each data type:
        # Impedance:         ZXX, ZXY, ZYX, ZYY
        # Impedance_Tipper:  ZXX, ZXY, ZYX, ZYY, TZX, TZY
        # Rho_Phs:           RhoXX, PhsXX, RhoXY, PhsXY, RhoYX, PhsYX, RhoYY, PhsYY
        # Rho_Phs_Tipper:    RhoXX, PhsXX, RhoXY, PhsXY, RhoYX, PhsYX, RhoYY, PhsYY,
        #                    RealTZX, ImagTZX, RealTZY, ImagTZY
        elseif contains(cline, "DataComp")
            tmp = split(cline)
            nDc = parse(Int, tmp[end])
            dataComp = Array{String}(nDc)
            nd = 0
            while nd < nDc
                cline = strip(readline(fid))
                nd = nd + 1
                dataComp[nd] = cline
            end

        # data values and standard errors
        elseif contains(cline, "Data Block")
            tmp   = split(cline)
            nData = parse(Int, tmp[end])
            nd    = 0
            rxID  = zeros(Int, nData)
            dcID  = zeros(Int, nData)
            freqID = zeros(Int, nData)
            if isComplex
                obsData = zeros(Complex128, nData)
            else
                obsData = zeros(Float64, nData)
            end
            dataErr = zeros(nData)

            while nd < nData
                cline = strip(readline(fid))
                while isempty(cline) || cline[1] == '#'
                    cline = strip(readline(fid))
                end
                cline = split(cline)
                nd = nd + 1
                freqID[nd] = parse(Int, cline[1])
                rxID[nd]   = parse(Int, cline[2])
                dcID[nd]   = parse(Int, cline[3])
                if isComplex
                    obsData[nd] = float(cline[4]) + float(cline[5]) * 1im
                    dataErr[nd] = float(cline[6])
                else
                    obsData[nd] = float(cline[4])
                    dataErr[nd] = float(cline[5])
                end

            end

        end

    end

    close(fid)

    datInfo = DataInfo(rxLoc, freqs, dataType, dataComp, rxID, freqID, dcID)

    return datInfo, obsData, dataErr

end
