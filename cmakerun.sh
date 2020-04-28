#!/bin/sh

# for armeabi-v7a and arm64-v8a cross compiling
toolchain_file="/media/shan/OS/Hiruna/temp-old/android-sdk-linux/ndk-bundle/build/cmake/android.toolchain.cmake"

# to create a artic library
sed -i 's/int main(int argc/int init_artic(int argc/g' src/main.cpp || die "main not found"
#touch src/interface.h
#echo "int init_artic(int argc, char ** argv);" > src/interface.h

echo "Creating build directories"

mkdir -p lib
rm -rf lib
mkdir lib
mkdir lib/armeabi-V7a
mkdir lib/arm64-v8a
mkdir lib/x86

mkdir -p build
rm -rf build
mkdir build
cd build

#echo "Building for x86 ..."
#
# for architecture x86
#cmake .. -DDEPLOY_PLATFORM=x86
#make -j 8

#rm -r ./*

#echo "Building for armeabi-V7a ..."
#
## for architecture armeabi-V7a
#cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE:STRING=$toolchain_file -DANDROID_PLATFORM=android-21 -DDEPLOY_PLATFORM:STRING="armeabi-v7a" -DANDROID_ABI="armeabi-v7a"
#ninja
#
#cp libartic_c.so ../lib/armeabi-V7a
#
#rm -r ./*
#
#echo "Building for arm64-v8a ..."
#
## for architecture arm66-v8a
cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE:STRING=$toolchain_file -DANDROID_PLATFORM=android-23 -DDEPLOY_PLATFORM:STRING="arm64-v8a" -DANDROID_ABI="arm64-v8a"
ninja
#
#cp libartic_c.so ../lib/arm64-v8a
#
#echo "\nShared libraries can be found in the lib directory"