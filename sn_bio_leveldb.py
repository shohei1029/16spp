#!/usr/bin/env python3

import leveldb
import pickle

# Created on 2016.12.24 by Shohei Nagata.
# leveldbの技術検証用スクリプト

db = leveldb.LevelDB('../leveldb')

test = "world"
in_sd = pickle.dumps(test)
print(in_sd)
db.Put('hello'.encode('utf-8'), in_sd)

out_sd = db.Get('hello')
out = pickle.loads(out_sd)
print(out)
