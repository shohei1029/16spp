#!/usr/bin/env python3

import sys

import toml

argvs = sys.argv

if __name__ == "__main__":
    with open(argvs[1]) as tomlfh:
        toml_d = toml.loads(tomlfh.read())

    for k in toml_d.keys():
        sys.stdout.write("{} ".format(k))
    
