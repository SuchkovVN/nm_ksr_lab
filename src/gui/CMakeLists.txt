set(TARGET_NAME nm_gui)

file(GLOB_RECURSE SOURCES "*.cc" "*.h")

add_executable(${TARGET_NAME} ${SOURCES}
    src/mainwindow.ui src/helpform.ui)

target_link_libraries(${TARGET_NAME} PRIVATE
    Qt6::Core 
    Qt6::PrintSupport
    Qt6::Widgets
    nm::core
    qcustomplot
)

set_target_properties(${TARGET_NAME} PROPERTIES
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

target_include_directories(${TARGET_NAME} PRIVATE include)

set(IMAGE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src/image.png" CACHE FILEPATH "path to help image")

target_compile_definitions(${TARGET_NAME} PUBLIC IMAGE_PATH="${IMAGE_PATH}")



