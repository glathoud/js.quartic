#!/usr/bin/env bash

index.scm
echo '_upload.sh: uploading "quartic" article...'
find . -name '*~' -exec rm {} \;
echo "mkdir js.quartic
cd js.quartic
mput -rf *.scm *.html *.js *.TXT *.xcf *.jpg
" | ncftp glat
echo '_upload.sh: upload successful!'
