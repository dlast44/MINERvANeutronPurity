add_library(obj NeutCands.cpp)
target_link_libraries(obj ${ROOT_LIBRARIES})
install(TARGETS obj DESTINATION lib)
install(FILES NeutCands.h DESTINATION include)
