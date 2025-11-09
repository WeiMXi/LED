message(STATUS "Looking for GLoBES")

set(LED_GLOBES_MINIMUM_REQUIRED 3.0.0)

if(NOT LED_BUILTIN_GLOBES)
    find_package(GLoBES ${LED_GLOBES_MINIMUM_REQUIRED})
    if(NOT GLobes_FOUND)
        set(LED_BUILTIN_GLOBES ON)
        message(NOTICE "***Notice: GLoBES not found (minimum required is ${LED_GLOBES_MINIMUM_REQUIRED}). For the time turning on LED_BUILTIN_GLOBES")
    endif()
endif()

if(LED_BUILTIN_GLOBES)
    message(STATUS "LED will use built-in GLoBES")
    # check built-in version
    if(LED_BUILTIN_GLOBES_VERSION VERSION_LESS LED_GLOBES_MINIMUM_REQUIRED)
        message(NOTICE "***Notice: Provided LED_BUILTIN_GLOBES_VERSION is ${LED_BUILTIN_GLOBES_VERSION}, which is less than the requirement (${LED_GLOBES_MINIMUM_REQUIRED}). Changing to ${LED_GLOBES_MINIMUM_REQUIRED}")
        set(LED_BUILTIN_GLOBES_VERSION ${LED_GLOBES_MINIMUM_REQUIRED})
    endif()
    # set download dest and URL
    set(LED_BUILTIN_GLOBES_SRC_DIR "${LED_PROJECT_3RDPARTY_DIR}/GLoBES-main")
    set(LED_BUILTIN_GLOBES_URL1 "https://github.com/WeiMXi/globes-cmake/archive/refs/heads/main.zip")
    set(LED_BUILTIN_GLOBES_URL2 "https://gh-proxy.com/https://github.com/WeiMXi/globes-cmake/archive/refs/heads/main.zip")
    # reuse or download
    include(FetchContent)
    if(EXISTS "${LED_BUILTIN_GLOBES_SRC_DIR}/CMakeLists.txt")
        FetchContent_Declare(GLoBES SOURCE_DIR "${LED_BUILTIN_GLOBES_SRC_DIR}")
        message(STATUS "Reusing GLoBES source ${LED_BUILTIN_GLOBES_SRC_DIR}")
    else()
        FetchContent_Declare(GLoBES SOURCE_DIR "${LED_BUILTIN_GLOBES_SRC_DIR}"
                                     URL "${LED_BUILTIN_GLOBES_URL1}" "${LED_BUILTIN_GLOBES_URL2}")
        message(STATUS "GLoBES will be downloaded from ${LED_BUILTIN_GLOBES_URL1} (with fallback to GitHub proxy) to ${LED_BUILTIN_GLOBES_SRC_DIR}")
    endif()
    # configure it
    message(STATUS "Downloading (if required) and configuring GLoBES (version: ${LED_BUILTIN_GLOBES_VERSION})")
    FetchContent_MakeAvailable(GLoBES)
    message(STATUS "Downloading (if required) and configuring GLoBES (version: ${LED_BUILTIN_GLOBES_VERSION}) - done")
endif()

if(NOT LED_BUILTIN_GLOBES)
    message(STATUS "Looking for GLoBES - found (version: ${GLOBES_VERSION})")
else()
    message(STATUS "Looking for GLoBES - built-in (version: ${LED_BUILTIN_GLOBES_VERSION})")
endif()
