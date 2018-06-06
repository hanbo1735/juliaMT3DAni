workspace()

push!(LOAD_PATH,"D:\\HB_Projects\\realEM\\MT3D_ANI_GeoPaper")

using MT3DFileIO

# (c) HB, Oct., 2016

println("=== Testing functions in module MT3DFileIO ===")

#datafile = "pek1dmodela_31Freqs.dat"
datafile = "100sec.dat"
modfile  = "pek1dmodela.mod"

#
println("Testing function readMT3DData ...")
@time datInfo, obsData, dataErr = readMT3DData(datafile)

#
println("Testing function readModel3DAni ...")
@time emMesh =  readModel3DAni(modfile)

#
println("Testing function writeMT3DData ...")
wdatfile = "wData.dat"
writeMT3DData(wdatfile, datInfo, obsData, dataErr)

println("Testing function writeModel3DAni ...")
wmodfile = "wModel.mod"
writeModel3DAni(wmodfile, emMesh)


#

println("=== All functions in module MT3DFileIO passed the test === ")
