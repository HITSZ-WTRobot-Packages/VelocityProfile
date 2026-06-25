add_library(VelocityProfileSCurve STATIC
    "./src/s_curve.cpp"
)

target_include_directories(VelocityProfileSCurve
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_options(VelocityProfileSCurve PRIVATE
    -Ofast
    -ffast-math
)

target_compile_definitions(VelocityProfileSCurve PUBLIC
    ARM_MATH_CM4
)

# link dependencies if any
target_link_libraries(VelocityProfileSCurve PUBLIC VelocityProfile::Core)

# alias for external use
add_library(VelocityProfile::SCurve ALIAS VelocityProfileSCurve)
