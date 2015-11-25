#!/bin/bash
find . -name "*.asciidoc" -exec asciidoc -a icons '{}' \;
echo "
You have run this:
find . -name '*.asciidoc' -exec asciidoc -a icons '{}' \;
then after connecting to web.sf.net:
cd htdocs
put *html
exit
"

sftp alpapan,jamg@web.sf.net

