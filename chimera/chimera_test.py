import os
#import math
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

rc("reset")
rc("open /Users/NagataShohei/Documents/bio-study/16spp/data/pdb/1rev.pdb")
rc("background solid white")
rc("select protein")
rc("color red selected")
rc("turn z -30")
rc("copy file /Users/NagataShohei/Documents/bio-study/16spp/scripts/chimera/chimera_out.png dpi 300 supersample 3")
rc("close all")
