add_executable("LinearElastic_executable" "${CMAKE_CURRENT_LIST_DIR}/test/test.cpp")
target_link_libraries("LinearElastic_executable" ${PROJECT_NAME})
add_test(NAME "LinearElastic" COMMAND "LinearElastic_executable")
