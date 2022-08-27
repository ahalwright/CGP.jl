#! /bin/bash
# This bash shell script uses the Linix tool "sed" to replace a line of files that create CSV files
#   to a version that runs on both Linux and MacOS
sed -i'.bak' 's@chomp(open("/etc/hostname") do f read(f,String) end)@readchomp(`hostname`)@' $1
