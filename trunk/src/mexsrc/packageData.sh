#!/bin/bash

printenv | grep G4 | grep DATA | sed 's/G4\w*=//' | xargs -I@ cp -r @ /var/MATLAB/G4Data/