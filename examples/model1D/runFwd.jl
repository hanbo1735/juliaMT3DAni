workspace()

# you can set environment variables under Julia as follows:
# ENV["OMP_NUM_THREADS"] = 4
# ENV["MKL_NUM_THREADS"] = 4
# ENV["MKL_PARDISO_OOC_MAX_CORE_SIZE"] = 50000


using MT3DFwdSolver
using MT3DFileIO
using YeeGridOperators
using MUMPS


#----------------------- read data file and model file ------------------------#

# directory of the input files
# inFileDir = pwd()*"\\"        # for Windows
inFileDir = pwd()*"/"      # for Linux

datafile = inFileDir * "31Freqs_RhoPhs.dat"
println("Reading data file $(datafile) ...")
@time datInfo, obsData, dataErr = readMT3DData(datafile)

modelfile = inFileDir * "Pek1DModel_grid2.mod"
println("Reading model file $(modelfile) ...")
@time emMesh = readModel3DAni(modelfile)



#------------------------ perform forward computation -------------------------#
# by defult, it will compute the according to input data information
# @time fwdData, eField, fwdDecomp = MT3DFwdComp(emMesh, datInfo, linSolver="mumps")

# you can specify the argument 'task' to perform a pure forward computation
@time fwdResp, = MT3DFwdComp(emMesh, datInfo, linSolver="mumps", task="pf")

# options for 'linSolver': mumps, pardiso, qmr, bicgstb. (case insensitive)


#------------------------------ output results --------------------------------#
# directory of the output files
outFileDir = identity(inFileDir)

# write forward response file
respfile = outFileDir * "31Freqs_RhoPhs.resp"
println("Writing out the forward responses ...")
writeMT3DResp(respfile, datInfo, fwdResp)

# write forward data file
# outdatafile = outFileDir * "31Freqs_RhoPhs_out.dat"
# writeMT3DData(outdatafile, datInfo, fwdData)


# write forward field file
# fieldfile = "31Freqs_RhoPhs.field"
# println("Writing out the edge fields ...")
# writeEdgeFields(fieldfile, eField, emMesh.gridSize)



println("=== Finishing forward problem ===")
