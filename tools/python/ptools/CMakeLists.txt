set(cmake-command "
  execute_process(
    COMMAND ${PYTHON_INTERPRETER} -m pip install --prefix=$ENV{DESTDIR}/${CMAKE_INSTALL_PREFIX} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
")

install(CODE "${cmake-command}")
