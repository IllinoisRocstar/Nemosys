# Setting executable files #####################################################

set(CARTESIAN2DMESH_SRCS cartesian2DMesh.C)

# Building executable ##########################################################

add_executable(cartesian2DMesh ${CARTESIAN2DMESH_SRCS})
target_link_libraries(cartesian2DMesh PRIVATE
    meshTools
    triSurface
    finiteVolume
    edgeMesh
)
target_link_libraries(cartesian2DMesh PRIVATE cfmesh)
target_compile_definitions(cartesian2DMesh PRIVATE NoSizeType)
