target_sources(${PROJECT_NAME} PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ContDynamics.tpp
        ${CMAKE_CURRENT_LIST_DIR}/LinearSys.tpp
        ${CMAKE_CURRENT_LIST_DIR}/NonlinearSys.tpp
        ${CMAKE_CURRENT_LIST_DIR}/polyReachableSetElement.tpp
        ${CMAKE_CURRENT_LIST_DIR}/polyReachOptions.tpp
        ${CMAKE_CURRENT_LIST_DIR}/ReachableSet.tpp
        ${CMAKE_CURRENT_LIST_DIR}/ReachableSetElement.tpp
        ${CMAKE_CURRENT_LIST_DIR}/TimeRelated.tpp)