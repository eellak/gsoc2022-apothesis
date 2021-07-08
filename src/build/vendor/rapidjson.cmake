include(ExternalProject)
ExternalProject_Add(
rapidjson
PREFIX "vendor/rapidjson"
GIT_REPOSITORY "https://github.com/Tencent/rapidjson.git"
GIT_TAG b557259f8813bed2d79e83bd2f92eb673579a82c 
TIMEOUT 10
CMAKE_ARGS
	-DRAPIDJSON_BUILD_TESTS=OFF
	-DRAPIDJSON_BUILD_DOC=OFF
	-DRAPIDJSON_BUILD_EXAMPLES=OFF
CONFIGURE_COMMAND ""
BUILD_COMMAND ""
INSTALL_COMMAND ""
UPDATE_COMMAND ""
)

ExternalProject_Get_Property(rapidjson source_dir)
set(RAPIDJSON_INCLUDE_DIR ${source_dir}/include)
