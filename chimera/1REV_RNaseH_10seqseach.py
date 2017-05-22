import os
#import math
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

#テスト用
#結局GUIでやるようにした．

rc("reset")
rc("open /Users/NagataShohei/Documents/bio-study/16spp/data/pdb/1rev.pdb")
rc("background solid white")
rc("select protein")
rc("color red selected")
rc("turn z -30")
#rc("copy file /Users/NagataShohei/tmp/hoge.png dpi 300 supersample 3")
rc("copy file ./out_chimera.png dpi 500 supersample 3")
rc("close all")
