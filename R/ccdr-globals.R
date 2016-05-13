#
#  ccdr-globals.R
#  ccdri
#
#  Created by Bryon Aragam (local) on 2/4/15.
#  Copyright (c) 2014-2015 Bryon Aragam (local). All rights reserved.
#

#
# PACKAGE CCDR: Global variables
#
#   CONTENTS:
#     .MACHINE_EPS
#

.MACHINE_EPS <- .Machine$double.eps ^ 0.5 # approximately 1.5e-08, used to be 1e-12
