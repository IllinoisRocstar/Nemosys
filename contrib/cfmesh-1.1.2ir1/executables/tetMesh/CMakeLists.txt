# Setting executable files #####################################################

set(TETMESH_SRCS tetMesh.C)

# Building executable ##########################################################

add_executable(tetMesh ${TETMESH_SRCS})
target_link_libraries(tetMesh PRIVATE
    meshTools
    triSurface
    finiteVolume
    edgeMesh
)
target_link_libraries(tetMesh PRIVATE cfmesh)
target_compile_definitions(tetMesh PRIVATE NoSizeType)
