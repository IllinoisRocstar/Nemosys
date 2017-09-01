#!/bin/bash

for item in $(ls *.inp); do
	vol2planeTransfer $item
done
