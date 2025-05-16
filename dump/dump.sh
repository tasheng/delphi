#!/bin/bash
pgm=dump

export DELLIBS=`dellib skelana dstana pxdst vfclap vdclap ux tanagra ufield bsaurus herlib trigger uhlib dstana`
export CERNLIBS=`cernlib  genlib packlib kernlib ariadne herwig jetset74`

# run patchy
ycmd=`which nypatchy`
command="$ycmd - ${pgm}.f $pgm.cra $pgm.ylog - - - ${pgm}.f .go"
echo "Running $command"
eval $command

# compile
for ending in .f .F ; do
    ls *$ending >/dev/null 2>&1
    if [ $? -eq 0 ]; then
	for file in *$ending  ; do
	    $FCOMP $FFLAGS -c $file
	done
    fi
done


$FCOMP $LDFLAGS *.o -o $pgm.exe $ADDLIB $DELLIBS ${CERNLIBS// -lnsl/}

# execute
./$pgm.exe 1>$pgm.log 2>$pgm.err

# cleanup
#rm -f *.f *.c *.o
