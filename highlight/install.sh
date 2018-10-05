#!/bin/bash
source=$(cd $(dirname "$0"); pwd)/gtksourceview/playmol.lang

destination=();
while IFS= read -d $'\0' -r file ; do
#  echo $file
  if echo $file | grep -q usr; then
    file=$(echo $file | grep usr | sed 's/def.lang/playmol.lang/g')
#    echo $file
    destination=("${destination[@]}" "$file");
  fi
done < <(locate -0 'language-specs/def.lang')

n=${#destination[@]}

if [ -z "$destination" ]; then
    echo "GtkSourceView was not found"
    exit
fi

if (( n > 1 )); then
    echo
    echo "[!] More than one GtkSourceView candidates found:"
    echo
    i=0
    for option in ${destination[@]}; do
      echo $i '=>' $option;
      ((i++));
    done
    echo
    echo "[?] Input integer corresponding to desired installation path"
    echo
    read j
#    echo $j
else
    j=0
fi

command="$(which cp) -rf $source ${destination[$j]}"
echo $command
eval $command
