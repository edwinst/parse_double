@echo off

rem Compiler options used:
rem     /nologo ... do not display Microsoft startup banner
rem     /Dxxx ... define preprocessor macros
rem     /EHa-s-c- ... turn off all exception handling
rem     /MT ... link to static run-time library
rem     /O2 ... optimize for speed
rem     /W3 ... "production quality" warnings
rem     /Zi ... put debug info in .pdb (not in .obj file)
rem     /link ... pass following options to the linker
rem Linker options used:
rem     /DEBUG:FULL ... move all debug info into the .pdb
rem     /DEBUG:NONE ... don't generate debug info
rem     /OUT:filename ... specifies the output file to create
rem     /INCREMENTAL:NO ... do not link incrementally
rem     /SUBSYSTEM:CONSOLE ... mark this as a Windows console application

set CXX_FLAGS=/nologo /W3 /EHa-s-c- /MT

cl %CXX_FLAGS% /Zi parse_double.cpp histogram.cpp /link /DEBUG:FULL /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /OUT:parse_double_debug.exe

cl %CXX_FLAGS% /DNDEBUG /O2 parse_double.cpp histogram.cpp /link /DEBUG:NONE /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /OUT:parse_double.exe
