#################################################################
# Search C files recursively in all subdirs
#################################################################

FILE(GLOB_RECURSE LIB_SRC *.c)
SET(LIB_SRC ${LIB_SRC} PARENT_SCOPE)

#################################################################
# Library target
#################################################################
ADD_LIBRARY(${LIB} ${LIB_SRC})

TARGET_LINK_LIBRARIES(${LIB} ${${PROJECT_NAME}_EXTERNAL_PROJECTS})

SET_TARGET_PROPERTIES(${LIB} PROPERTIES VERSION ${${LIB}_VERSION} SOVERSION ${${LIB}_SOVERSION})
EXPORT(TARGETS ${LIB} APPEND FILE ${ROOT_OUTPUT_PATH}/${LIB}Targets.cmake)

