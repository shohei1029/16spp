#!/usr/bin/env python3

import sys

#memo
# Attrs の ID欄にドメインやセクションなどの名前をそれぞれ追加していく
# 2016.3.4 by Shohei N.

#hard coding
#variations = ["-int_sec", "-RT_sec", "-RNase_sec", "-prot_sec"] #including separators
#variations = ["-Protease", "-p51_RT", "-p15", "-Integrase"] #including separators
variations = ["-RVP", "-RVT_1", "-RVT_thumb", "-RVT_connect", "-RNase_H", "-Integrase_Zn", "-rve", "-IN_DBD_C"] #including separators

if __name__ == '__main__':
    for line in sys.stdin:
        header =  line
        sys.stdout.write(header)
        break
    
    for line in sys.stdin:
        id =  line.split("\t")[0]
        rest =  line.split("\t")[1:]
    
        for s in variations:
            ls = []
            id2 = id + s
            ls.append(id2)
            ls.extend(rest)
            sys.stdout.write("\t".join(ls))


