#!/bin/sh
chmod +x ./amd64-linux/php-cgi
exec ./amd64-linux/php-cgi -c ./amd64-linux/php-cgi.ini "$@"