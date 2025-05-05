#!/bin/bash


#export generator=zgpy
#export generator=wwpy
#export generator="."
#export generator=eegg
#export generator=kora
#export generator=wphact24_cc
#export generator=wphact21_nc
#export generator=kk2f4144_qqar
#export generator=kk2f4143_qq
export generator=kk2f4144_ttsp
#export generator=qedbk23eegg

export year=1998
export version=v2

if [[ $year -eq 1998 ]]; then
    export tag="v98e1"
fi

if [[ $year -eq 1994 ]]; then
    export tag="v94c2"
fi

if [[ "$generator" != "." ]]; then
    folder="simulation"
else
    folder="collision_data"
fi


cd /eos/user/z/zhangj/DELPHI/${folder}/${year}_${version}/${tag}/${generator}/

ls ${generator}*.root > allfiles.txt

# 2. Split that list into chunks of 50 lines each
split -l 50 allfiles.txt chunk_

# 3. Loop over each chunk and hadd its files
i=1
for f in chunk_*; do
  out=$(printf "merged_%02d.root" "$i")
  echo "Merging $(wc -l < "$f") files into $out"
  hadd -f "$out" $(<"$f")
  i=$((i+1))
done
