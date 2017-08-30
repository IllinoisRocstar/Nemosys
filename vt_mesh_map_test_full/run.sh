#!/bin/bash

for item in $(ls *.inp); do
	vt_mesh_map $item
done
