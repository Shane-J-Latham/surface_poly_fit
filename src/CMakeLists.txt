IF (DEFINED VCPKG_INSTALLED_DIR)
  # List `vcpkg` dependencies from VCPKG_INSTALLED_DIR so we can install them
  # together with the python API.
  FILE(
    GLOB VCPKG_DEPENDENCIES
    RELATIVE "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/"
    "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/bin/*"
  )
ENDIF()

IF ( CGAL_FOUND )
    SET(
      spf_cgal_SRCS
      polyhedral_surface.cpp
      polyhedral_surface_py.cpp
      spf_cgal.cpp
      spf_cgal_py.cpp
    )

    ADD_LIBRARY(_spf_cgal MODULE ${spf_cgal_SRCS})
    TARGET_LINK_LIBRARIES(
      _spf_cgal
      CGAL::CGAL
      ${PYTHON_LIBRARIES}
    )
    PYTHON_EXTENSION_MODULE(_spf_cgal)

    # Install the `vcpkg` dependencies that are used in our ${PROJECT_NAME}
    # library in the same DESTINATION as the python API. Note that it
    # only supports collecting the runtime dependencies for Windows, Linux and
    # macOS platforms.
    IF (DEFINED SKBUILD)
      INSTALL(
        TARGETS _spf_cgal
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_DEPENDENCIES}
        PRE_EXCLUDE_REGEXES ".*"
        DESTINATION surface_poly_fit
    )
    ELSE()
      MESSAGE(STATUS "PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR=${PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR}")
      INSTALL(
        TARGETS _spf_cgal
        RUNTIME_DEPENDENCIES
        PRE_INCLUDE_REGEXES ${VCPKG_DEPENDENCIES}
        PRE_EXCLUDE_REGEXES ".*"
        DESTINATION ${PCSR_PYTHON_SITE_PACKAGES_INSTALL_DIR}
      )
    ENDIF()
ENDIF()

