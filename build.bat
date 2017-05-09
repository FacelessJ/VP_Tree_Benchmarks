@echo OFF
if not exist "build" mkdir build

set buildType=Release
if "%~1"=="Debug" (
	set buildType=Debug
	)

cd build
echo Build configuration: %buildType%
cmake -DCMAKE_BUILD_TYPE=%buildType% ..
cmake --build .
cd ..