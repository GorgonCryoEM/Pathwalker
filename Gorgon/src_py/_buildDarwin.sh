#!/bin/bash

cp ../../../libpyGORGON.dylib libpyGORGON.dylib
cp ../../../source/ExternalLibraries/boost_1_35_0/bin.v2/libs/python/build/darwin/release/macosx-version-10.4/threading-multi/libboost_python-mt-1_35.dylib libboost_python-mt-1_35.dylib 

rm libpyGORGON.so
rm -R dist
python _buildExecutableDarwin.py py2app
cp libpyGORGON.dylib dist/gorgon.app/Contents/Frameworks
cp libboost_python-mt-1_35.dylib dist/gorgon.app/Contents/Resources/ 


cd dist/gorgon.app/Contents/Resources/lib/python2.5/lib-dynload/
ln -s ../../../../Frameworks/libpyGORGON.dylib libpyGORGON.so
chmod +rx libpyGORGON.so
ln -s ../OpenGL OpenGL
chmod +rx OpenGL
cd ../../../../../../../


ln -s libpyGORGON.dylib libpyGORGON.so

