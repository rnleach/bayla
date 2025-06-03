#!/bin/bash -l

PROJDIR="$(dirname $(realpath $0))"
LIBSRCDIR="$PROJDIR/lib"

cd $PROJDIR

echo
echo "------------------------------------------------------------------------------------------------------------------------"
echo "                                            ********** Build Test **********"
echo

CC=cc

CFLAGS="-Wall -Werror -Wno-unknown-pragmas -march=native -fPIC -std=c11"
CFLAGS="$CFLAGS -D_DEFAULT_SOURCE"
CFLAGS="$CFLAGS -I$LIBSRCDIR"

LDLIBS="-lm"

if [ "$#" -gt 0 -a "$1" = "debug" ]
then
    echo "debug build"
    CFLAGS="$CFLAGS -O0 -DCOY_PROFILE -D_MAG_TRACK_MEM_USAGE"
    CFLAGS="$CFLAGS -g -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function"
elif [ "$#" -gt 0 -a "$1" != "clean" -o \( "$#" = 0 \) ]
then
    echo "release build"
    CFLAGS="$CFLAGS -O3 -DNDEBUG"
fi

if [ "$#" -gt 0 -a "$1" = "clean" ] 
then
    echo "clean compiled test program"
    echo
    rm -rf bayla-test *.dSYM
fi

if [ "$#" -gt 0 -a "$1" = "debug" ]
then
    echo
    echo "build bayla-test"
    $CC $CFLAGS test/main.c -o bayla-test $LDLIBS

    if [ "$#" -gt 1 -a "$2" = 'test' ]
    then
        echo "run bayla-test" && echo && echo && ./bayla-test
    fi
fi

if [ "$#" -gt 0 -a "$1" = "test" ]
then
    echo
    echo "build bayla-test"
    $CC $CFLAGS test/main.c -o bayla-test $LDLIBS && echo && echo "run bayla-test" && echo && echo && ./bayla-test
fi

echo
echo "                                              ********** Done Test **********"
echo "------------------------------------------------------------------------------------------------------------------------"
echo

