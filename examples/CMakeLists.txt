# Source files in this directory
set(SRC_FILES
    5thGKdV
    7thorder
    Eckhaus
    example1
    example2
    KdVB
    kerrNLS
    KS
    NLS
)

# Add executables and link to library
if(GINACDE_BUILD_SHARED_LIBS OR GINACDE_BUILD_STATIC_LIBS)
    foreach(src ${SRC_FILES})
        add_executable(${src} ${src}.cpp)
        target_link_libraries(${src}
            PUBLIC $<IF:$<BOOL:GINACDE_BUILD_SHARED_LIBS>,GiNaCDE,GiNaCDE_static>
        )
    endforeach()
endif()

# Add meta-target to compile
add_custom_target(examples)
add_dependencies(examples ${SRC_FILES})
