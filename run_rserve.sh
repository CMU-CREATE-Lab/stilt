#!/bin/bash

cd /usr/local/stilt
/usr/bin/R CMD Rserve --vanilla --RS-conf Rserve.conf >> rserve.log 2>&1
