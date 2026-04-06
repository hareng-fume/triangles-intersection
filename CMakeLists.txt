cmake_minimum_required(VERSION 3.23)

project(triangle-proj LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSTIONS OFF)

add_executable(triangle-app
	triangle_intersect.cpp
	triangle_intersect.h
	main.cpp
)

if(MSVC)
	target_compile_options(triangle-app PRIVATE /W4 /WX)
else()
	target_compile_options(triangle-app PRIVATE -Wall -Werror)
endif()

