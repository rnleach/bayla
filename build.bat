rem
rem build.bat clean - deletes test file
rem build.bat - builds the compiler "optimized" version of the test and runs it.
rem build.bat debug - builds the "debug" version of the test and does NOT launch it.
rem

@echo off
cls
SETLOCAL

rem
rem Set up compiler flags
rem
SET common_flags=/std:c11 /TC /utf-8 /nologo /favor:INTEL64 /arch:AVX2 /W3 /WX
IF "%1"=="debug" (GOTO Debug) ELSE (GOTO Release)

:Debug
@echo Debug Build
SET flags=/Od /Zi /DCOY_PROFILE
GOTO Operation

:Release
@echo Release Build
SET flags=/O2 /DNDEBUG /DCOY_PROFILE
GOTO Operation

rem
rem The main operation....build or clean?
rem
:Operation
IF "%1"=="clean" GOTO Clean 
GOTO BuildAll

rem
rem Clean up operations
rem
:Clean
@echo Clean
del *.exe *.obj *.pdb *.ilk *.dll
GOTO EndSuccess

rem
rem Build All
rem
:BuildAll

@echo Build
cl %common_flags% %flags% test\main.c /link /OUT:bayla-test.exe
IF "%1"=="test" (GOTO Test) ELSE (GOTO EndSuccess)

rem
rem Test
rem
:Test
@echo(
@echo ---------- Running Test ----------
@echo(
bayla-test.exe
@echo(
@echo ---------- Finished Test ----------
@echo(


rem
rem Exit Quietly with no errors
rem
:EndSuccess

