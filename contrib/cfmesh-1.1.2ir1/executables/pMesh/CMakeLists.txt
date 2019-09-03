# Setting executable files #####################################################

set(PMESH_SRCS pMesh.C)

# Building executable ##########################################################

add_executable(pMesh ${PMESH_SRCS})
target_link_libraries(pMesh PRIVATE
    meshTools
    triSurface
    finiteVolume
    edgeMesh
)
target_link_libraries(pMesh PRIVATE cfmesh)
target_compile_definitions(pMesh PRIVATE NoSizeType)
