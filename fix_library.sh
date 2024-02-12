#!/bin/zsh

install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libccp4c.8.0.dylib $CCP4/lib/libccp4c.8.0.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libclipper-ccp4.2.dylib $CCP4/lib/libclipper-ccp4.2.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libclipper-contrib.2.dylib $CCP4/lib/libclipper-contrib.2.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libclipper-core.2.dylib $CCP4/lib/libclipper-core.2.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libclipper-minimol.2.dylib $CCP4/lib/libclipper-minimol.2.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libclipper-mmdb.2.dylib $CCP4/lib/libclipper-mmdb.2.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libfftw.2.dylib $CCP4/lib/libfftw.2.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/libmmdb2.0.dylib $CCP4/lib/libmmdb2.0.dylib sails
install_name_tool -change /Users/jenkins/workspace/CCP4/series-8.0-osx/devtools/install/lib/librfftw.2.dylib $CCP4/lib/librfftw.2.dylib sails
