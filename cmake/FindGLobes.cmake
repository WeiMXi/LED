message(STATUS "Looking for GLoBES")

set(LIV_GLOBES_MINIMUM_REQUIRED 3.0.0)

if(NOT LIV_BUILTIN_GLOBES)
    find_package(GLoBES ${LIV_GLOBES_MINIMUM_REQUIRED})
    if(NOT GLobes_FOUND)
        set(LIV_BUILTIN_GLOBES ON)
        message(NOTICE "***Notice: GLoBES not found (minimum required is ${LIV_GLOBES_MINIMUM_REQUIRED}). For the time turning on LIV_BUILTIN_GLOBES")
    endif()
endif()

if(LIV_BUILTIN_GLOBES)
    message(STATUS "LIV will use built-in GLoBES")
    # check built-in version
    if(LIV_BUILTIN_GLOBES_VERSION VERSION_LESS LIV_GLOBES_MINIMUM_REQUIRED)
        message(NOTICE "***Notice: Provided LIV_BUILTIN_GLOBES_VERSION is ${LIV_BUILTIN_GLOBES_VERSION}, which is less than the requirement (${LIV_GLOBES_MINIMUM_REQUIRED}). Changing to ${LIV_GLOBES_MINIMUM_REQUIRED}")
        set(LIV_BUILTIN_GLOBES_VERSION ${LIV_GLOBES_MINIMUM_REQUIRED})
    endif()
    # set download dest and URL
    set(LIV_BUILTIN_GLOBES_SRC_DIR "${LIV_PROJECT_3RDPARTY_DIR}/GLoBES-main")
    set(LIV_BUILTIN_GLOBES_URL1 "https://github.com/WeiMXi/globes-cmake/archive/refs/heads/main.zip")
    set(LIV_BUILTIN_GLOBES_URL2 "https://gh-proxy.com/https://github.com/WeiMXi/globes-cmake/archive/refs/heads/main.zip")
    # reuse or download
    include(FetchContent)
    if(EXISTS "${LIV_BUILTIN_GLOBES_SRC_DIR}/CMakeLists.txt")
        FetchContent_Declare(GLoBES SOURCE_DIR "${LIV_BUILTIN_GLOBES_SRC_DIR}")
        message(STATUS "Reusing GLoBES source ${LIV_BUILTIN_GLOBES_SRC_DIR}")
    else()
        FetchContent_Declare(GLoBES SOURCE_DIR "${LIV_BUILTIN_GLOBES_SRC_DIR}"
                                     URL "${LIV_BUILTIN_GLOBES_URL1}" "${LIV_BUILTIN_GLOBES_URL2}")
        message(STATUS "GLoBES will be downloaded from ${LIV_BUILTIN_GLOBES_URL1} (with fallback to GitHub proxy) to ${LIV_BUILTIN_GLOBES_SRC_DIR}")
    endif()
    # configure it
    message(STATUS "Downloading (if required) and configuring GLoBES (version: ${LIV_BUILTIN_GLOBES_VERSION})")
    FetchContent_MakeAvailable(GLoBES)
    message(STATUS "Downloading (if required) and configuring GLoBES (version: ${LIV_BUILTIN_GLOBES_VERSION}) - done")
endif()

if(NOT LIV_BUILTIN_GLOBES)
    message(STATUS "Looking for GLoBES - found (version: ${GLOBES_VERSION})")
else()
    message(STATUS "Looking for GLoBES - built-in (version: ${LIV_BUILTIN_GLOBES_VERSION})")
endif()
