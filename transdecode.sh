#!/bin/bash

t=$1

~/applications/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $t
~/applications/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t $t
