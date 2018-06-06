workspace()

# ENV["OMP_NUM_THREADS"] = 1
# ENV["MKL_NUM_THREADS"] = 1

# ENV["MKL_PARDISO_OOC_MAX_CORE_SIZE"] = 50000

#push!(LOAD_PATH, "/home/hanbo/code/MT3D_ANI_GeoPaper")
push!(LOAD_PATH,"D:\\code\\MT3D_ANI_GeoPaper")


using MT3DFwdSolver
using MT3DFileIO
using YeeGridOperators

using MUMPS
# using Pardiso

# data file
datafile = "pek1dmodela_7Freqs.dat"
println("Reading data file $(datafile) ...")
@time datInfo, obsData, dataErr = readMT3DData(datafile)

# model file
modelfile  = "pek1dmodela.mod"
println("Reading model file $(modelfile) ...")
@time emMesh = readModel3DAni(modelfile)


@time fwdResp, = MT3DFwdComp(emMesh, datInfo, linSolver="qmr", task="pf", iterTol=1e-6)
respfile = "pek1dmodela_7Freqs.resp"
println("Writing out the forward responses ...")
writeMT3DResp(respfile, datInfo, fwdResp)

#=
fieldfile = "field.field"
println("Writing out the edge fields ...")
writeEdgeFields(fieldfile, eField, emMesh.gridSize)
=#

#=
@time predData, = MT3DFwdComp(emMesh, datInfo)
pdfile = "pek1dmodela_pred.dat"
writeMT3DData(pdfile, datInfo, predData)
=#

println("=== Finishing forward problem ===")
