# Setting executable files #####################################################

set(CARTESIANMESH_SRCS cartesianMesh.C)

# Building executable ##########################################################

add_executable(cartesianMesh ${CARTESIANMESH_SRCS})
target_link_libraries(cartesianMesh PRIVATE
    meshTools
    triSurface
    finiteVolume
    edgeMesh
)
target_link_libraries(cartesianMesh PRIVATE cfmesh)
target_compile_definitions(cartesianMesh PRIVATE NoSizeType)
