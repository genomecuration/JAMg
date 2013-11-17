#!/bin/bash
rm index.asciidoc~ -f
asciidoc --backend slidy index.asciidoc
firefox index.html
