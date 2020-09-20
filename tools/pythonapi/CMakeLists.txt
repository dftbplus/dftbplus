set(cmake-command "
  execute_process(
    COMMAND ${PYTHON_INTERPRETER} setup.py install --prefix=$ENV{DESTDIR}/${CMAKE_INSTALL_PREFIX}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
")

install(CODE "${cmake-command}")
