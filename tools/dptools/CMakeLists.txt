set(prefix ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX})
set(cmake-command "
  execute_process(
    COMMAND ${PYTHON_INTERPRETER} -m pip install --no-deps --prefix ${prefix} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
")

install(CODE "${cmake-command}")
