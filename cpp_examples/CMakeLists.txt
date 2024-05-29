function(add_cam_example ccfile)
  get_filename_component(file ${ccfile} NAME_WE)
  add_executable(${file} ${ccfile})
  target_link_libraries(${file} PUBLIC CAM)
  target_include_directories(${file} PUBLIC ${PROJECT_SOURCE_DIR})
  add_test(NAME "${file}_cxx" COMMAND "${file}" WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endfunction()

add_cam_example(cpp_example.cxx)
