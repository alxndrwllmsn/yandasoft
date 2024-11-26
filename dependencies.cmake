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
set ( LOFAR_COMMON_TAG             f35117f80d12c9bd4a363c23fb00c8e0c49b1cee      CACHE  STRING    "lofar-common version"             FORCE )
set ( LOFAR_BLOB_TAG               tags/1.3.1      CACHE  STRING    "lofar-blob version"               FORCE )
set ( BASE_ASKAP_TAG               ffd9370a6e5ebc5e3f07e2383bca2df431dd8269    CACHE  STRING    "base-askap version"               FORCE )
set ( BASE_IMAGEMATH_TAG           9ce81c9021c013eb1dbdce83fe6e94a6228d68b9    CACHE  STRING    "base-imagemath version"           FORCE )
set ( BASE_ASKAPPARALLEL_TAG       1ccc9663cf58a71007db5dc48e4ad3beb0a4d419    CACHE  STRING    "base-askapparallel version"       FORCE )
set ( BASE_SCIMATH_TAG             feb949247c760972ae3c7a0936bb7869ee3ebc35    CACHE  STRING    "base-scimath version"             FORCE )
set ( BASE_ACCESSORS_TAG           f4fd6bb7ac565f82d97caf715c7f97b34854bbf5    CACHE  STRING    "base-accessors version"           FORCE )
set ( BASE_COMPONENTS_TAG          ffcdab6a4de43812b9beb0dc90cdd5e2f7502bab    CACHE  STRING    "base-components version"          FORCE )
# TOS related repos are not versioned yet, so pinned this build with commit hash values.
