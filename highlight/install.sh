#!/bin/bash
source=$(cd $(dirname "$0"); pwd)/playmol.lang
destination=$(locate language-specs/def.lang | sed 's/def.lang/playmol.lang/g')
if [ -z "$destination" ]; then
    echo "GtkSourceView was not found"
    exit
fi
command="$(which cp) -rf $source $destination"
echo $command
eval $command

