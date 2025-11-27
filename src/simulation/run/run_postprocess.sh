#!/bin/bash

set -e


postProcess -func forces

postProcess -func "pressureCoefficient"

foamToVtk -lastTime