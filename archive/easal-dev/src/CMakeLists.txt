target_sources(easal
    PRIVATE
          	${CMAKE_CURRENT_LIST_DIR}/main.cpp
    PUBLIC
)

if(WITH_QT)
    ADD_DEFINITIONS(-DQT_NO_KEYWORDS)
    add_library(easalqtgui "")
    add_subdirectory(gui_Qt)
    target_link_libraries(easalqtgui LINK_PUBLIC easalcore Qt5::Widgets)
    if(WITH_CAF)
            target_link_libraries(easalqtgui
                    LINK_PRIVATE
                    ${CAF_LIBRARIES})
    endif()


endif()


add_library(easalcore "")
add_subdirectory(easalcore)

if(WITH_CAF)
    target_link_libraries(easalcore
       LINK_PRIVATE
       ${CAF_LIBRARIES})
endif()


