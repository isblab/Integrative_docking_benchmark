#include(GoogleTest)
find_package(GTest REQUIRED) 

add_executable(tests main.cpp
                     src/PathFinderTest.cpp
		     #src/ActiveConstraintRegionTest.cpp
		     #src/AtlasNodeTest.cpp
		     #src/CayleyPointTest.cpp
                     #src/sanity_tests.cpp
                     #src/jacobian_tests.cpp 
                     #src/feature_tests.cpp 
                     )

include_directories(../src)

target_link_libraries(tests 
                GTest::GTest easalcore 
	    	${GFLAGS}/lib/libgflags${CMAKE_SHARED_LIBRARY_SUFFIX}
                ${GLOG}/lib/libglog${CMAKE_SHARED_LIBRARY_SUFFIX})
gtest_discover_tests(tests TEST_PREFIX Testing)
