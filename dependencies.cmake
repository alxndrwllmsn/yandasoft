# A list of specific version for the dependencies of this build.
#
# The contents here will aid in 'pinning' the contents of the build. Changing this file constitutes changing the build
# and this means that the integration version should be promoted. The version numbers (tags) specified here will be
# injected into the CMakeLists.txt file in the External_project_add definitions and will, therefore, control the
# versions of various sub libs built and integrated. They need not be official versions and can be arbitrary branches or
# tags such as dev feature branches etc.
#
# You can affect the build content by manipulating the "tag/branch" column. Should be something that git expects for that 
# repo.
#     Identifier                   tag/branch       cache  type      description                        force it
set ( ASKAP_CMAKE_TAG              tags/1.5.0      CACHE  STRING    "askap.cmake tools"                FORCE )
set ( ASKAP_DEV_TAG                tags/2.13.0     CACHE  STRING    "askap dev tools"                  FORCE )
set ( LOFAR_COMMON_TAG             ae919181bb2adbcf1141ea01f64340ea36c83685    CACHE  STRING    "lofar-common version"             FORCE )
set ( LOFAR_BLOB_TAG               tags/1.3.1      CACHE  STRING    "lofar-blob version"               FORCE )
set ( BASE_ASKAP_TAG               044f3ef17c4209e39d45a32763c3edc22563dbb4    CACHE  STRING    "base-askap version"               FORCE )
set ( BASE_LOGFILTERS_TAG          ca5f91a431a0df668199a9c6f05383bc0c38e135    CACHE  STRING    "base-logfilters version"          FORCE )
set ( BASE_IMAGEMATH_TAG           df34e16364d5cb77334ae25cfcf21ecf93e8dac4    CACHE  STRING    "base-imagemath version"           FORCE )
set ( BASE_ASKAPPARALLEL_TAG       7a08893da050e4c710984eaa93cc1f9208eca5c2    CACHE  STRING    "base-askapparallel version"       FORCE )
set ( BASE_SCIMATH_TAG             12e2cd852e58d23740dd8b06b472ed6d82a3ddb6    CACHE  STRING    "base-scimath version"             FORCE )
set ( BASE_ACCESSORS_TAG           016599d64669fd8ad5b71d272b2b2f07f7ccb6d6    CACHE  STRING    "base-accessors version"           FORCE )
set ( BASE_COMPONENTS_TAG          96024d1b9cdce2e6e29685c09c60c0c28f67d5b5    CACHE  STRING    "base-components version"          FORCE )
