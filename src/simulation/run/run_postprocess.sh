#!/bin/bash

set -ex


postProcess -func forces

postProcess -func "pressureCoefficient"

foamToVtk -lastTime